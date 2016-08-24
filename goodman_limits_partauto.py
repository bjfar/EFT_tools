""" Translating limits from Goodman et. al. (2011) 
    i.e. Majorana fermion dark matter
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy import interp

import old.limit_comparison.translate_couplings as t
import old.limit_comparison.comparison_tools as tools
import old.limit_comparison.parameters as p

rootpath = "/home/farmer/mathematica/DMFormFactor_13086288/"

#=====================================

# coupling for operator M1,M4
# g_{\chi-q}
# I assume they mean they coupling is set equal for all quarks
# (or rather M is the same for all quarks, then scaled by each of their masses for some FCNC reason)
# note, for M2 and M3, g = i*g1q
def g1q(mq,M):
  return mq / (2*M**3)
g2q = g1q  # can ignore i while not worrying too much about interference
g3q = g1q
g4q = g1q

#coupling for operator M5,M6
def g5q(M):
   return 1 / (2*M**2)
g6q = g5q

# coupling for operator M7,M9
# Note, for M8 and M10, g = i*g7G
# g_{\chi-G}
def g7G(M):
  return p.alpha_s / (8*M**3)
g8G  = g7G # again ignoring i's for now
g9G  = g7G
g10G = g7G

# We have limits on M as function of m_\chi
# Need to translate to limits on quark/gluon couplings, and then nucleon couplings
# Goodman et. al. say "We will assume that the interaction is dominated by
# only one of the above operators in the table", which I will take to mean that
# for the M1 limit, ALL the quarks couplings are contributing, but not the gluon
# couplings, and vice versa for M7.

#================================
# Goodman limit curves
M1_mchi, M1_limit = np.loadtxt("goodman_limits_data/M1+3_tevatron.csv",delimiter=",").T
M2_mchi, M2_limit = np.loadtxt("goodman_limits_data/M2+4_tevatron.csv",delimiter=",").T
M5_mchi, M5_limit = np.loadtxt("goodman_limits_data/M5+6_tevatron.csv",delimiter=",").T
M7_mchi, M7_limit = np.loadtxt("goodman_limits_data/M7+9_tevatron.csv",delimiter=",").T
M8_mchi, M8_limit = np.loadtxt("goodman_limits_data/M8+10_tevatron.csv",delimiter=",").T

# Do some (constant) extrapolation on the low end, and interpolation in the middle.
smooth_M1_mchi = np.logspace(0,np.log10(np.max(M1_mchi)),num=50)
smooth_M2_mchi = np.logspace(0,np.log10(np.max(M2_mchi)),num=50)
smooth_M5_mchi = np.logspace(0,np.log10(np.max(M5_mchi)),num=50)
smooth_M7_mchi = np.logspace(0,np.log10(np.max(M7_mchi)),num=50)
smooth_M8_mchi = np.logspace(0,np.log10(np.max(M8_mchi)),num=50)
smooth_M1_limit = interp(smooth_M1_mchi, M1_mchi, M1_limit)
smooth_M2_limit = interp(smooth_M2_mchi, M2_mchi, M2_limit)
smooth_M5_limit = interp(smooth_M5_mchi, M5_mchi, M5_limit)
smooth_M7_limit = interp(smooth_M7_mchi, M7_mchi, M7_limit)
smooth_M8_limit = interp(smooth_M8_mchi, M8_mchi, M8_limit)

# Visually check limit digitization and interpolation
fig = plt.figure(figsize=(5,3))
ax = fig.add_subplot(111)
#ax.plot(M1_mchi, M1_limit, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
#ax.plot(M2_mchi, M2_limit, label="$\overline{\chi}\gamma_5\chi\overline{q}q$ (Tevatron)", lw=2)
#ax.plot(M5_mchi, M5_limit, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma^\mu q$ (Tevatron)", lw=2)
#ax.plot(M7_mchi, M7_limit, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
#ax.plot(M8_mchi, M8_limit, label="$\overline{\chi}\gamma_5\chi G G$ (Tevatron)", lw=2)
#ax.plot(smooth_M1_mchi, smooth_M1_limit, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
#ax.plot(smooth_M2_mchi, smooth_M2_limit, label="$\overline{\chi}\gamma_5\chi\overline{q}q$ (Tevatron)", lw=2)
#ax.plot(smooth_M5_mchi, smooth_M5_limit, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma^\mu q$ (Tevatron)", lw=2)
#ax.plot(smooth_M7_mchi, smooth_M7_limit, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
#ax.plot(smooth_M8_mchi, smooth_M8_limit, label="$\overline{\chi}\gamma_5\chi G G$ (Tevatron)", lw=2)
# Alternate Labels
ax.plot(smooth_M1_mchi, smooth_M1_limit, label="$M1,M3$", lw=2)
ax.plot(smooth_M2_mchi, smooth_M2_limit, label="$M2,M4$", lw=2)
ax.plot(smooth_M5_mchi, smooth_M5_limit, label="$M5,M6$", lw=2)
ax.plot(smooth_M7_mchi, smooth_M7_limit, label="$M7,M9$", lw=2)
ax.plot(smooth_M8_mchi, smooth_M8_limit, label="$M8,M10$", lw=2)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$ (GeV)")
ax.set_ylabel(r"$M_{*}$ (GeV)")
plt.legend(frameon=False)
plt.tight_layout()
fig.savefig("Mstar_lim_check.png")

# If it looks ok, replace original data with interpolated/extrapolated data
M1_mchi  = smooth_M1_mchi 
M2_mchi  = smooth_M2_mchi 
M5_mchi  = smooth_M5_mchi 
M7_mchi  = smooth_M7_mchi 
M8_mchi  = smooth_M8_mchi 
M1_limit = smooth_M1_limit
M2_limit = smooth_M2_limit
M5_limit = smooth_M5_limit
M7_limit = smooth_M7_limit
M8_limit = smooth_M8_limit

# Fill in degenerate limits
M3_mchi, M3_limit = M1_mchi, M1_limit
M4_mchi, M4_limit = M2_mchi, M2_limit
M6_mchi, M6_limit = M5_mchi, M5_limit
M9_mchi, M9_limit = M7_mchi, M7_limit
M10_mchi, M10_limit = M8_mchi, M8_limit

#=========================================

# Manual matching of quark/gluon couplings to nucleon couplings
# Store limit curves in list of lists of dictionaries which stores all info needed
# to generate plots

Nops = 15 + 1 # plus one to offset the zero indexing
limits_O = [[] for i in range(Nops)]

#=========================================
# Limits on O1 from M1 and M7 (tevatron) 
# i.e. scalar coupling to quarks and gluons
# ---> scalar coupling to nucleons

g1u_limit = g1q(p.mu,M1_limit)
g1d_limit = g1q(p.md,M1_limit)
g1s_limit = g1q(p.ms,M1_limit)
g1c_limit = g1q(p.mc,M1_limit)
g1b_limit = g1q(p.mb,M1_limit)
g1t_limit = g1q(p.mt,M1_limit)

g7G_limit = g7G(M7_limit)

# limits on nucleon couplings from quark couplings (M1):
g1p_M1q = t.calc_cN_scalar(g1u_limit,g1d_limit,g1s_limit,g1c_limit,g1b_limit,g1t_limit,0.,p.data_p)
g1n_M1q = t.calc_cN_scalar(g1u_limit,g1d_limit,g1s_limit,g1c_limit,g1b_limit,g1t_limit,0.,p.data_n)

# limits on nucleon couplings from gluon couplings (M7):
g1p_M7G = t.calc_cN_scalar(0.,0.,0.,0.,0.,0.,g7G_limit,p.data_p)
g1n_M7G = t.calc_cN_scalar(0.,0.,0.,0.,0.,0.,g7G_limit,p.data_n)

limits_O[1]+= [{"data" : g1n_M1q**2 * p.mWeak**4, 
                "mchi" : M1_mchi, 
                "label": "$\overline{\chi}\chi\overline{q}q$ (Tevatron)"}
              ,{"data" : g1n_M7G**2 * p.mWeak**4, 
                "mchi" : M7_mchi, 
                "label": "$\overline{\chi}\chi G G$ (Tevatron)"}]

#=========================================
# Limits on O4 from M6 (tevatron) 
# i.e. axial-vector coupling to quarks
# ---> axial-vector coupling to nucleons

g6_limit = g6q(M6_limit) # same for u,d,s

# limits on nucleon couplings from quark couplings (M6):
# again assume same coupling for all quarks
g4p_M6q = t.calc_cN_axialvector(g6_limit,g6_limit,g6_limit,p.data_p)
g4n_M6q = t.calc_cN_axialvector(g6_limit,g6_limit,g6_limit,p.data_n)

limits_O[4]+= [{"data" : g4n_M6q**2 * p.mWeak**4, 
                "mchi" : M6_mchi, 
                "label": "$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)"}]

#=========================================
# Limits on O6 from M4 and M10 (tevatron) 
# i.e. pseudoscalar coupling to quarks and gluons
# ---> pseudoscalar coupling to nucleons
# (with coupling to pseudoscalar WIMP current)

g4u_limit = g4q(p.mu,M4_limit)
g4d_limit = g4q(p.md,M4_limit)
g4s_limit = g4q(p.ms,M4_limit)
g4c_limit = g4q(p.mc,M4_limit)
g4b_limit = g4q(p.mb,M4_limit)
g4t_limit = g4q(p.mt,M4_limit)

g10G_limit = g10G(M10_limit)

# limits on nucleon couplings from quark couplings (M1):
# (extra factor is from mapping of O_4^R --> O_6^NR
g6p_M4q = (-p.mN/M4_mchi) * t.calc_cN_pseudoscalar(g4u_limit,g4d_limit,g4s_limit,g4c_limit,g4b_limit,g4t_limit,0.,p.data_p)
g6n_M4q = (-p.mN/M4_mchi) * t.calc_cN_pseudoscalar(g4u_limit,g4d_limit,g4s_limit,g4c_limit,g4b_limit,g4t_limit,0.,p.data_n)

# limits on nucleon couplings from gluon couplings (M7):
g6p_M10G = (-p.mN/M10_mchi) * t.calc_cN_pseudoscalar(0.,0.,0.,0.,0.,0.,g10G_limit,p.data_p)
g6n_M10G = (-p.mN/M10_mchi) * t.calc_cN_pseudoscalar(0.,0.,0.,0.,0.,0.,g10G_limit,p.data_n)

limits_O[6]+= [{"data" : g6n_M4q**2 * p.mWeak**4, 
                "mchi" : M4_mchi, 
                "label": r"$\overline{\chi}\gamma_5\chi\overline{q}\gamma_5q$ (Tevatron)"}
              ,{"data" : g6n_M10G**2 * p.mWeak**4, 
                "mchi" : M10_mchi, 
                "label": r"$\overline{\chi}\gamma_5\chi G \tilde{G}$ (Tevatron)"}]

#==================================================
# Limits on O10 from M3 and M9  (tevatron) 
# i.e. pseudoscalar coupling to quarks
# ---> pseudoscalar coupling to nucleons
# (with coupling to scalar WIMP current)

g3u_limit = g3q(p.mu,M3_limit)
g3d_limit = g3q(p.md,M3_limit)
g3s_limit = g3q(p.ms,M3_limit)
g3c_limit = g3q(p.mc,M3_limit)
g3b_limit = g3q(p.mb,M3_limit)
g3t_limit = g3q(p.mt,M3_limit)

g9G_limit = g9G(M9_limit)

# Not sure this really makes sense, should probably set all couplings
# equal to largest(smallest?) of limits. Though actually I think this is assuming
# the same M*, i.e. mediator mass / coupling, for each effecting operator, so
# maybe ok.
g10p_M3q = t.calc_cN_pseudoscalar(g3u_limit,g3d_limit,g3s_limit,g3c_limit,g3b_limit,g3t_limit,0,p.data_p)
g10n_M3q = t.calc_cN_pseudoscalar(g3u_limit,g3d_limit,g3s_limit,g3c_limit,g3b_limit,g3t_limit,0,p.data_n)

g10p_M9G = t.calc_cN_pseudoscalar(0,0,0,0,0,0,g9G_limit,p.data_p)
g10n_M9G = t.calc_cN_pseudoscalar(0,0,0,0,0,0,g9G_limit,p.data_n)

limits_O[10]+= [{"data": g10n_M3q**2 * p.mWeak**4, 
                 "mchi": M3_mchi, 
                 "label": r"$\overline{\chi}\chi\overline{q}\gamma_5 q$ (Tevatron)"}
               ,{"data": g10n_M9G**2 * p.mWeak**4, 
                 "mchi": M9_mchi,  
                 "label": r"$\overline{\chi}\chi G \tilde{G}$ (Tevatron)"}]

#=========================================

# Automated extraction of estimated Xenon100 limits in all EFT couplings
excluded = [0,2]
for i in range(Nops):
   if i not in excluded:
      print "=== Extracting estimated Xenon100 225 live day limit for operator {0} ===".format(i)
      Xenon100_est_EFT_lim   = tools.get_c_curve("c{0}p=c{0}n".format(i)) #equal proton/neutron couplings
      Xenon100_est_EFT_lim_n = tools.get_c_curve("c{0}n".format(i))       #neutron-only couplings
      limits_O[i] += [{"data" : np.array(Xenon100_est_EFT_lim)**2,  
                       "mchi" : tools.masses, 
                       "label": "Xenon100 N={0} EFT (c{1}p=c{1}n)".format(tools.Nevents,i)}
                     ,{"data" : np.array(Xenon100_est_EFT_lim_n)**2,
                       "mchi" : tools.masses, 
                       "label": "Xenon100 N={0} EFT (c{1}n)".format(tools.Nevents,i)}]

# Automated plotting of various limits in all EFT couplings
def add_curve(ax,datadict):   
   ax.plot( datadict["mchi"]
          , datadict["data"]
          , label=datadict["label"]
          , lw=2) 

Nplots = np.sum([len(limits_O[i])>0 for i in range(Nops)])
fig = plt.figure(figsize=(6,Nplots*4))
nextfree = 1
for i in range(Nops):
   if len(limits_O[i])>0:
      ax = fig.add_subplot(Nplots,1,nextfree)
      nextfree+=1
      for curve in limits_O[i]:
         add_curve(ax,curve)
      ax.set_xscale("log")
      ax.set_yscale("log")
      ax.set_xlabel("$m_\chi$")
      ax.set_ylabel(r"$g_{{{0}}}^2 \times\, m^4_\mathrm{{weak}}$".format(i))
      plt.legend(frameon=False)

plt.tight_layout()
fig.savefig("tevatron_lims_auto.png")

