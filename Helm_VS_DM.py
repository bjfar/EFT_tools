# Comparison of recoil spectra using Helm form factors against the density matrix treatment

import numpy as np
import matplotlib.pyplot as plt
import ROOT
import old.limit_comparison.translate_couplings as tc

fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)

#helm = np.loadtxt("recoil_spectrum_tables/Xe131/Xe131_NR_HelmFF_c1p=c1n.dat").T
dens = np.loadtxt("recoil_spectrum_tables/Xe131/Xe131_NR_c1p=c1n.dat").T

print dens.shape

masses = np.array([3,4,5,6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,100,200,300,500,700,800,1000,2000])
#masses = np.array([5,10,50])
col = np.where(masses==10)[0][0] + 1 # column of the mass we want to extract
print col

# Spectra from Ludwig
f = ROOT.TFile("run10_10_nom_Ludwig.root", 'r')
#print f.ls() # see what is in the file
g = f.Get("gEnr10")
ludx = np.fromiter(g.GetX(), dtype='d', count=g.GetN())
ludy = np.fromiter(g.GetY(), dtype='d', count=g.GetN())

# Spectra from Boris
#f1 = ROOT.TFile("dr_de_Boris.root",'r')
#print f1.ls() # see what is in the file
#g1 = f1.Get("dRdE_Mass10_Leff-1.0000_Qy-1.0000_LCE1.0000NR_c1_isoscalar")
#oldx = np.fromiter(g1.GetX(), dtype='d', count=g1.GetN())
#oldy = np.fromiter(g1.GetY(), dtype='d', count=g1.GetN())

# Scale Ludwig's spectra to match the Helm FF version, at some arbitrary point
h10i = (dens[0] == 1)
l10i = ( np.abs(ludx-1) < 5e-3 )
print l10i
K =  ludy[l10i] / dens[col][h10i]
print K

# # Scale Ludwig's spectra according to theoretical correspondence.
# # Should be 10^-44 cm^2 SI cross-section
# lud_ref = 1e-44
# 
# # First need raw dimensionful coupling, i.e. without extra weak-scale factors
# c1_ref = 1 # reference coupling normalisation used to generate dRdE
# mWeak = 246.2 # GeV
# mWIMP = 10 # GeV
# g1sq = c1_ref**2 * mWeak**-4
# sigmaSI = g1sq * tc.gsq_to_sigma(mWIMP,1)
# print "sigmaSI for c=1:", sigmaSI
# 
# # Ok need to scale our results so that match 10^-44 sigmaSI. Should go as ratio of couplings squared, i.e. cross-section ratio.
# K = lud_ref / sigmaSI
# print K


# # Validation of conversion against CDMSII Si limits, EFT vs normal, for 10 GeV WIMP, SI/O1
# CDMSIISi_sigmaSI_10GeV = 2.17e-41 # 1304.3706v3 figure 4, also abstract
# CDMSIISi_c1sqXmWeak4 = 3e-3 # 1503.03379 table 1
# 
# g1sq_CDMSIISi = CDMSIISi_c1sqXmWeak4 * mWeak**-4
# CDMSIISi_sigmaSI_EFT = g1sq_CDMSIISi * tc.gsq_to_sigma(10,1)
# 
# print "CDMSIISi sigmaSI 10 GeV limit, from standard analysis VS converted from EFT analysis"
# print CDMSIISi_sigmaSI_10GeV
# print CDMSIISi_sigmaSI_EFT

#ax.plot(helm[0], helm[col]*K, c="r",label="Helm FF")
ax.plot(dens[0], dens[col]*K, c="g",label="Density M.")
ax.plot(ludx,ludy,c="b",label="Ludwig")
#ax.set_xscale("log")
ax.set_yscale("log")
#ax.set_ylim(1,1e8)
ax.set_xlim(0,16)
ax.set_ylabel("dR/dE")
ax.set_xlabel("Recoil energy (keV)")
ax.legend(loc=1, frameon=False, framealpha=0,prop={'size':10})
fig.savefig("Helm_VS_DM_10GeV_O1.png")
