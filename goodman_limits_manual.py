""" Translating limits from Goodman et. al. (2011) 
    i.e. Majorana fermion dark matter
"""

import numpy as np
import matplotlib.pyplot as plt
#from scipy.interpolate import interp1d
from scipy import interp

import limit_comparison.translate_couplings as t
import limit_comparison.comparison_tools as tools
import limit_comparison.parameters as p

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

#coupling for operator M5,M6
def g5q(M):
   return 1 / (2*M**2)
g6q = g5q

# coupling for operator M7,M9
# Note, for M8 and M10, g = i*g7G
# g_{\chi-G}
def g7G(M):
  return p.alpha_s / (8*M**3)
g9G = g7G

# We have limits on M as function of m_\chi
# Need to translate to limits on quark/gluon couplings, and then nucleon couplings
# Goodman et. al. say "We will assume that the interaction is dominated by
# only one of the above operators in the table", which I will take to mean that
# for the M1 limit, ALL the quarks couplings are contributing, but not the gluon
# couplings, and vice versa for M7.

#================================
# Goodman limit curves
M1_mchi, M1_limit = np.loadtxt("goodman_limits/M1+3_tevatron.csv",delimiter=",").T
M6_mchi, M6_limit = np.loadtxt("goodman_limits/M5+6_tevatron.csv",delimiter=",").T
M7_mchi, M7_limit = np.loadtxt("goodman_limits/M7+9_tevatron.csv",delimiter=",").T

# Do some (constant) extrapolation on the low end, and interpolation in the middle.
smooth_M1_mchi = np.logspace(0,np.log10(np.max(M1_mchi)),num=50)
smooth_M6_mchi = np.logspace(0,np.log10(np.max(M6_mchi)),num=50)
smooth_M7_mchi = np.logspace(0,np.log10(np.max(M7_mchi)),num=50)
smooth_M1_limit = interp(smooth_M1_mchi, M1_mchi, M1_limit)
smooth_M6_limit = interp(smooth_M6_mchi, M6_mchi, M6_limit)
smooth_M7_limit = interp(smooth_M7_mchi, M7_mchi, M7_limit)

# Visually check limit digitization and interpolation
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(M1_mchi, M1_limit, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
ax.plot(M6_mchi, M6_limit, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)", lw=2)
ax.plot(M7_mchi, M7_limit, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
ax.plot(smooth_M1_mchi, smooth_M1_limit, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
ax.plot(smooth_M6_mchi, smooth_M6_limit, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)", lw=2)
ax.plot(smooth_M7_mchi, smooth_M7_limit, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$M_{*}$")
plt.legend(frameon=False)
plt.tight_layout()
fig.savefig("MstarM1M6M7_lims.png")

# If it looks ok, replace original data with interpolated/extrapolated data
M1_mchi  = smooth_M1_mchi 
M6_mchi  = smooth_M6_mchi 
M7_mchi  = smooth_M7_mchi 
M1_limit = smooth_M1_limit
M6_limit = smooth_M6_limit
M7_limit = smooth_M7_limit

# Fill in degenerate limits
M3_mchi, M3_limit = M1_mchi, M1_limit
M9_mchi, M9_limit = M7_mchi, M7_limit

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

# Translate to limits on nucleon couplings
#-------

# limits on nucleon couplings from quark couplings (M1):
g1p_M1q = t.calc_cN_scalar(g1u_limit,g1d_limit,g1s_limit,g1c_limit,g1b_limit,g1t_limit,0.,p.data_p)
g1n_M1q = t.calc_cN_scalar(g1u_limit,g1d_limit,g1s_limit,g1c_limit,g1b_limit,g1t_limit,0.,p.data_n)

# limits on nucleon couplings from gluon couplings (M7):
g1p_M7G = t.calc_cN_scalar(0.,0.,0.,0.,0.,0.,g7G_limit,p.data_p)
g1n_M7G = t.calc_cN_scalar(0.,0.,0.,0.,0.,0.,g7G_limit,p.data_n)

# Compare to cross-section computed using Goodman formulae (eq. 4,5)
sigmaSI_M1q_v2 = 4*t.Mred(M1_mchi,p.mN)**2 / np.pi * (0.082) * (1/(2*M1_limit**3))**2 * p.GeV2_to_cm2
sigmaSI_M7G_v2 = 4*t.Mred(M7_mchi,p.mN)**2 / np.pi * (5.0)   * (1/(8*M7_limit**3))**2 * p.GeV2_to_cm2

#=========================================
# Limits on O4 from M6 (tevatron) 
# i.e. axial-vector coupling to quarks
# ---> axial-vector coupling to nucleons

g6_limit = g6q(M6_limit) # same for u,d,s

# Translate to limits on nucleon couplings
#-------

# Summary of what is going on:
# Coupling plots:
#   (1) Get limit on M* based on M6 operator in Goodman et al.
#       - Convert it to limit on coupling of M6 operator (trivial)
#   (2) Get limit on sigma_SD^n from Xenon225 SD paper
#       - Convert it to limit on M* by inverting Goodman et al eq. 3
#       - ... and then to limit on coupling as in (1)
#   (3) Get limit on coupling from EFT N=5 predictions for Xenon225
# 
# Nucleon cross-section plots
#   (1) Convert limit on M* based on M6 to limit on sigma_SD^n using 
#       Goodman eq. 3
#   (2) Use limit on sigma_SD^n from Xenon225 paper directly
#   (3) Convert limit on coupling from EFT N=5 prediction into limit on 
#       sigma_SD^n using Goodman eq. 3

# limits on nucleon couplings from quark couplings (M6):
# again assume same coupling for all quarks
g4p_M6q = -4*t.calc_cN_axialvector(g6_limit,g6_limit,g6_limit,p.data_p)
g4n_M6q = -4*t.calc_cN_axialvector(g6_limit,g6_limit,g6_limit,p.data_n)

# Goodman formula for cross-section (eq. 5 in their longer paper)
# (includes conversion to cm^2)
sigmaSD_M6q_v2 = 9.18e-40 * t.Mred(M6_mchi,p.mN)**2 * (300./M6_limit)**4
# EDIT: Seems close, might be small (in log units) factor of maybe 4 off

# Note, limits from experiment assume coupling to *only* proton or neutron
# How to arrange this fundamental couplings?
# note: c_N = cu*Delu + cd*Deld + cs*Dels
# For proton:  Delp_u =  0.77, Delp_d = -0.4,  Delp_s = -0.12
# For neutron: Deln_u = -0.4,  Deln_d =  0.77, Delp_s = -0.12
# Pretty weird arrangement needed to cancel
# But, Xenon100 is much more sensitive to the neutron coupling anyway,
# so perhaps just compare the "equal coupling to all quarks" model against
# the constraints on neutron coupling.

#==================================================
# Limits on O10 from M3 and M9  (tevatron) 
# i.e. pseudoscalar coupling to quarks
# ---> pseudoscalar coupling to nucleons

g3u_limit = g3q(p.mu,M3_limit)
g3d_limit = g3q(p.md,M3_limit)
g3s_limit = g3q(p.ms,M3_limit)
g3c_limit = g3q(p.mc,M3_limit)
g3b_limit = g3q(p.mb,M3_limit)
g3t_limit = g3q(p.mt,M3_limit)

g9G_limit = g9G(M9_limit)

# Translate to limits on nucleon couplings
#-------

# Not sure this really makes sense, should probably set all couplings
# equal to largest(smallest?) of limits. Though actually I think this is assuming
# the same M*, i.e. mediator mass / coupling, for each effecting operator, so
# maybe ok.
g10p_M3q = t.calc_cN_pseudoscalar(g3u_limit,g3d_limit,g3s_limit,g3c_limit,g3b_limit,g3t_limit,0,p.data_p)
g10n_M3q = t.calc_cN_pseudoscalar(g3u_limit,g3d_limit,g3s_limit,g3c_limit,g3b_limit,g3t_limit,0,p.data_n)

g10p_M9G = t.calc_cN_pseudoscalar(0,0,0,0,0,0,g9G_limit,p.data_p)
g10n_M9G = t.calc_cN_pseudoscalar(0,0,0,0,0,0,g9G_limit,p.data_n)

#=========================================

# Get SI Xenon100 225 live day limits (1207.5988)
mX225, sigmaSI_X225 = np.loadtxt(rootpath+"/goodman_limits/xenon225livedays.csv",delimiter=",").T
sigmaSI_X225 = 10**sigmaSI_X225 # data in file is in log10 units

# Compare to SI Xenon100 estimated limits via EFT:
c1_lim_Xenon100 = tools.get_c_curve("c1p=c1n") # SI operator, equal proton/neutron couplings
c1n_lim_Xenon100 = tools.get_c_curve("c1n") # SI operator, neutron only couplings

# Convert to/from sigma_SD (in cm^-2)

# WIMP-nucleon couplings squared
g1nsq_M1q    = g1n_M1q**2
g1nsq_M7G    = g1n_M7G**2
g1sq_limEFT  = np.array(c1_lim_Xenon100)**2  * p.mWeak**-4
g1nsq_limEFT = np.array(c1n_lim_Xenon100)**2 * p.mWeak**-4
g1sq_X225    = sigmaSI_X225 / t.gsq_to_sigma(mX225,1)

# cross sections
sigmaSI_M1q     = g1nsq_M1q    * t.gsq_to_sigma(M1_mchi,1)
sigmaSI_M7G     = g1nsq_M7G    * t.gsq_to_sigma(M7_mchi,1)
sigmaSI_limEFT  = g1sq_limEFT  * t.gsq_to_sigma(tools.masses,1)
sigmaSIn_limEFT = g1nsq_limEFT * t.gsq_to_sigma(tools.masses,1)
#sigmaSI_X225   #nothing to do

N = 3
fig = plt.figure(figsize=(12,N*4))
ax = fig.add_subplot(N,2,1)
ax.plot(M1_mchi, sigmaSI_M1q,    label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2) 
ax.plot(M7_mchi, sigmaSI_M7G,    label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2) 
ax.plot(M1_mchi, sigmaSI_M1q_v2, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron) v.2", lw=2, c='purple') 
ax.plot(M7_mchi, sigmaSI_M7G_v2, label="$\overline{\chi}\chi G G$ (Tevatron) v.2", lw=2, c='m') 
ax.plot(tools.masses,  sigmaSI_limEFT, label="Xenon100 N=5 EFT (c1p=c1n)", c='r', lw=2)
ax.plot(tools.masses,  sigmaSIn_limEFT, label="Xenon100 N=5 EFT (c1n)", c='orange', lw=2)
ax.plot(mX225,   sigmaSI_X225,   label="Xenon100 225 live days", c='k', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$\sigma^N_{SI}$")
#ax.set_ylabel(r"$g_1^2 * (m^4_\mathrm{weak})$")
plt.legend(frameon=False)

# WIMP-nucleon cross-sections
ylims = ax.get_ylim() # Get y limits from cross-section plot
ax = fig.add_subplot(N,2,2)
# Note, plotting coupling to neutrons, but proton coupling is very similar unless quark/gluon-level couplings highly tuned
ax.plot(M1_mchi, g1nsq_M1q   * p.mWeak**4, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
ax.plot(M1_mchi, g1nsq_M7G   * p.mWeak**4, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
ax.plot(tools.masses,  g1sq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT (c1p=c1n)", c='r', lw=2)
ax.plot(tools.masses,  g1nsq_limEFT* p.mWeak**4, label="Xenon100 N=5 EFT (c1n)", c='orange', lw=2)
ax.plot(mX225,   g1sq_X225   * p.mWeak**4, label="Xenon100 225 live days", c='k', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$g_1^2 \times\, m^4_\mathrm{weak}$")
plt.legend(frameon=False)
# Add second axis for nucleon cross-section units
tools.add_second_scale(ax, t.gsq_to_sigmaOnMr2[1]*p.mWeak**-4, "$\sigma^N_{SI} / \mu_\chi^2$", ylims)

# PLOT FOR NOTE:
#-------
fig2 = plt.figure(figsize=(6,4))
ax = fig2.add_subplot(111)
ax.plot(M1_mchi, g1nsq_M1q   * p.mWeak**4, label="$\overline{\chi}\chi\overline{q}q$ (Tevatron)", lw=2)
ax.plot(M1_mchi, g1nsq_M7G   * p.mWeak**4, label="$\overline{\chi}\chi G G$ (Tevatron)", lw=2)
ax.plot(tools.masses,  g1sq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT (c1p=c1n)", c='r', lw=2)
ax.plot(tools.masses,  g1nsq_limEFT* p.mWeak**4, label="Xenon100 N=5 EFT (c1n)", c='orange', lw=2)
#ax.plot(mX225,   g1sq_X225   * p.mWeak**4, label="Xenon100 225 live days", c='k', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$c_1^2 \times\, m^4_\mathrm{weak}$")
plt.legend(frameon=False)
plt.tight_layout()
fig2.savefig("note_M1M7_lims.png")
#-------

#-- Now the SD plots ---

# Get SD Xenon100 225 live day limits (1301.6620)
mX225_SDn, sigmaSDn_X225 = np.loadtxt(rootpath+"/goodman_limits/xenon225livedays_SDn.csv",delimiter=",").T
#mX225_SDp, sigmaSDp_X225 = np.loadtxt(rootpath+"/goodman_limits/xenon225livedays_SDp.csv",delimiter=",").T
sigmaSDn_X225 = 10**sigmaSDn_X225 # data in file is in log10 units
#sigmaSDp_X225 = 10**sigmaSDn_X225 # data in file is in log10 units

# SD Xenon100 estimated limits via EFT:
c4_lim_Xenon100 = tools.get_c_curve("c4p=c4n") # SD operator, equal proton/neutron couplings
#c4p_lim_Xenon100 = tools.get_c_curve("c4p")    # SD operator, coupling only to proton
c4n_lim_Xenon100 = tools.get_c_curve("c4n")    # SD operator, coupling only to neutron; Will be basically same as c4p=c4n case due to small sensitivity of Xenon100 to c4p.

# Convert to/from sigma_SD (in cm^-2)

# WIMP-nucleon couplings squared
g4nsq_M6q    = g4n_M6q**2
g4sq_limEFT  = np.array(c4_lim_Xenon100)**2  * p.mWeak**-4
g4nsq_limEFT = np.array(c4n_lim_Xenon100)**2 * p.mWeak**-4
g4nsq_X225   = sigmaSDn_X225 / t.gsq_to_sigma(mX225_SDn,4)

# cross sections
sigmaSDn_M6q    = g4nsq_M6q    * t.gsq_to_sigma(M6_mchi,4)
sigmaSD_limEFT  = g4sq_limEFT  * t.gsq_to_sigma(tools.masses,4)
sigmaSDn_limEFT = g4nsq_limEFT * t.gsq_to_sigma(tools.masses,4)
#sigmaSDn_X225   #nothing to do


ax = fig.add_subplot(N,2,3)
ax.plot(M6_mchi, sigmaSDn_M6q,   label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)", lw=2, c='purple') 
ax.plot(M6_mchi, sigmaSD_M6q_v2, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron) v.2", lw=2, c='magenta') 
ax.plot(tools.masses,  sigmaSD_limEFT,  label="Xenon100 N=5 EFT (c4n=c4p)", c='r', lw=2)
ax.plot(tools.masses,  sigmaSDn_limEFT, label="Xenon100 N=5 EFT (c4n)", c='orange', lw=2)
ax.plot(mX225_SDn, sigmaSDn_X225, label="Xenon100 225 live days (n)", c='k', lw=2)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$\sigma^N_{SD}$")
#ax.set_ylabel(r"$c_1^2 * (m^4_\mathrm{weak})$")
plt.legend(frameon=False)

# Plot of coupling limits
ylims = ax.get_ylim() # Get y limits from cross-section plot
ax = fig.add_subplot(N,2,4)
ax.plot(M6_mchi, g4nsq_M6q    * p.mWeak**4, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)", lw=2, c='purple')
ax.plot(tools.masses,  g4sq_limEFT  * p.mWeak**4, label="Xenon100 N=5 EFT $(c_4^p=c_4^n)$", c='r', lw=2)
ax.plot(tools.masses,  g4nsq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT $(c_4^n)$", c='orange', lw=2)
ax.plot(mX225_SDn, g4nsq_X225 * p.mWeak**4, label="Xenon100 225 live days (n)", c='k', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$g_4^2 \times\, m^4_\mathrm{weak}$")
plt.legend(frameon=False)
# Add second axis for nucleon cross-section units
tools.add_second_scale(ax, t.gsq_to_sigmaOnMr2[4]*p.mWeak**-4, "$\sigma^N_{SD} / \mu_\chi^2$", ylims)

plt.tight_layout()

# PLOT FOR NOTE:
#-------
fig2 = plt.figure(figsize=(6,4))
ax = fig2.add_subplot(111)
ax.plot(M6_mchi, g4nsq_M6q    * p.mWeak**4, label="$\overline{\chi}\gamma_5\gamma_\mu\chi\overline{q}\gamma_5\gamma^\mu q$ (Tevatron)", lw=2, c='purple')
ax.plot(tools.masses,  g4sq_limEFT  * p.mWeak**4, label="Xenon100 N=5 EFT $(c_4^p=c_4^n)$", c='r', lw=2)
ax.plot(tools.masses,  g4nsq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT $(c_4^n)$", c='orange', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$c_4^2 \times\, m^4_\mathrm{weak}$")
ax.set_ylim(np.array(ax.get_ylim())/10.)
plt.legend(frameon=False)
plt.tight_layout()
fig2.savefig("note_M6_lims.png")
#-------


# Notes:
# Cross section and coupling plots don't quite seem consistent? Tevatron bounds appear stronger relative to the Xenon bounds in the cross-section plots, might still be some conversions not quite correct.

#--- Now try something new

# SD Xenon100 estimated limits via EFT:
c10_lim_Xenon100 = tools.get_c_curve("c10p=c10n") # SD operator, equal proton/neutron couplings
c10n_lim_Xenon100 = tools.get_c_curve("c10n")    # SD operator, coupling only to neutron

# Convert to/from sigma_O10 (in cm^-2)

# WIMP-nucleon couplings squared
g10psq_M3q   = g10p_M3q**2
g10nsq_M3q   = g10n_M3q**2

g10psq_M9G   = g10p_M9G**2
g10nsq_M9G   = g10n_M9G**2

g10sq_limEFT = np.array(c10_lim_Xenon100)**2  * p.mWeak**-4
g10nsq_limEFT= np.array(c10n_lim_Xenon100)**2 * p.mWeak**-4

# cross sections
# don't know how to do it yet
# sigmaSDn_M6q    = g4nsq_M6q    * t.gsq_to_sigma(M6_mchi,4)
# sigmaSD_limEFT  = g4sq_limEFT  * t.gsq_to_sigma(tools.masses,4)
# sigmaSDn_limEFT = g4nsq_limEFT * t.gsq_to_sigma(tools.masses,4)
# #sigmaSDn_X225   #nothing to do

# Plot of coupling limits
ylims = ax.get_ylim() # Get y limits from cross-section plot
ax = fig.add_subplot(N,2,6)
ax.plot(M3_mchi, g10nsq_M3q * p.mWeak**4, label="$\overline{\chi}\chi\overline{q}\gamma_5 q$ (Tevatron)", lw=2, c='b')
ax.plot(M9_mchi, g10nsq_M9G * p.mWeak**4, label="$\overline{\chi}\chi G \widetilde{G}$ (Tevatron)", lw=2, c='c')
ax.plot(tools.masses,  g10sq_limEFT  * p.mWeak**4, label="Xenon100 N=5 EFT $(c_{10}^p=c_{10}^n)$", c='r', lw=2)
ax.plot(tools.masses,  g10nsq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT $(c_{10}^n)$", c='orange', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$g_{10}^2 \times\, m^4_\mathrm{weak}$")
plt.legend(frameon=False)
# Add second axis for nucleon cross-section units
#tools.add_second_scale(ax, t.gsq_to_sigmaOnMr2[4]*p.mWeak**-4, "$\sigma^N_{SD} / \mu_\chi^2$", ylims)
plt.tight_layout()

# PLOT FOR NOTE:
#-------
fig2 = plt.figure(figsize=(6,4))
ax = fig2.add_subplot(111)
ax.plot(M3_mchi, g10nsq_M3q * p.mWeak**4, label="$\overline{\chi}\chi\overline{q}\gamma_5 q$ (Tevatron)", lw=2, c='b')
ax.plot(M9_mchi, g10nsq_M9G * p.mWeak**4, label="$\overline{\chi}\chi G \widetilde{G}$ (Tevatron)", lw=2, c='c')
ax.plot(tools.masses,  g10sq_limEFT  * p.mWeak**4, label="Xenon100 N=5 EFT $(c_{10}^p=c_{10}^n)$", c='r', lw=2)
ax.plot(tools.masses,  g10nsq_limEFT * p.mWeak**4, label="Xenon100 N=5 EFT $(c_{10}^n)$", c='orange', lw=2)
#ax.plot(mX225,   g1sq_X225   * p.mWeak**4, label="Xenon100 225 live days", c='k', lw=2)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("$m_\chi$")
ax.set_ylabel(r"$c_{10}^2 \times\, m^4_\mathrm{weak}$")
plt.legend(frameon=False)
plt.tight_layout()
fig2.savefig("note_M3M9_lims.png")
#-------



fig.savefig("tevatron_lims.png")

