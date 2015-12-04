import numpy as np
from parameters import GeV2_to_cm2, mN # nuclear data, conversion factors, etc.

import normalise_spectra as ns

rootpath = "/home/farmer/mathematica/DMFormFactor_13086288/"

# Extraction of estimated EFT couplings from mathematica-produced tables
#========================================
# Compare to Xenon100 estimated limits via EFT:
pathtodat = rootpath+"/EFTcoeffplotdata/"
Nevents = 5
exposure = 225*34 * 0.5  # 50% exposure to roughly account for cuts/acceptance
iso = "Comb" # use isotope abundance weighted ("Combined") recoil spectra
masses = np.array([5,6,10,30,50,100,300,500,700,800,1000,2000,3000,5000]) 

# Get estimated limit curve from EFT predictions
def get_c_curve(Op):
  c_lim_est = []
  for m in masses:
     print "Extracting: ",iso, Op, "mWIMP=",m
     data = np.loadtxt("{0}/Xe{1}/Xe{1}_{2}.dat".format(pathtodat,iso,Op))
     print "Normalising spectrum for Xe iso={0}, operator={1}, mWIMP={2}".format(iso,Op,m)
     index = np.where(masses==m)[0][0]+1 #first column is ER
     print index
     spectrum = data[:,[0,index]] 
     scaled_spectrum, coupling = ns.normalise_spectrum(spectrum,exposure,normrate=Nevents,normrange=(8,240),verbose=True)
     c_lim_est += [coupling]
  return c_lim_est  

# Plotting helpers
#==================================================
# Add second scale to axis
def add_second_scale(ax,K,label,ylims=None):
   ax2 = ax.twinx()
   ax2.set_xscale(ax.get_xscale())
   ax2.set_yscale(ax.get_yscale())
   ax2.set_ylabel(label)
   ax2.set_xlim(*ax.get_xlim())
   if ylims==None:
      y1, y2 = ax.get_ylim()
      ax2.set_ylim(y1*K, y2*K)
   else:
      # Automatically adjust the original axis limits to match the ones requested for the second axis
      ax2.set_ylim(ylims[0], ylims[1])
      ax.set_ylim(ylims[0]/K, ylims[1]/K)

