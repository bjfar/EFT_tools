import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import normalise_spectra as ns

pathtodat = "/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/"

labels = [["c{0}p".format(i),"c{0}n".format(i),"c{0}p=c{0}n".format(i)] for i in range(1,16)]

masses = np.array([6,10,50,100,500,1000])
isotopes = [128,129,130,131,132,134,136]

# Exposure to use for normalisation
exposure = 225*34

def get_norm_spectra(iso, lab, mass):
   print "Extracting: ",iso, lab, "mWIMP=",mass
   data = np.loadtxt("{0}/Xe{1}/Xe{1}_{2}.dat".format(pathtodat,iso,lab))
   print "Normalising spectrum for Xe iso={0}, operator={1}, mWIMP={2}".format(iso,lab,mass)
   index = np.where(masses==mass)[0][0]+1 #first column is ER
   spectrum = data[:,[0,index]]
   scaled_spectrum, c = ns.normalise_spectrum(spectrum,exposure,normrate=5)

   # output normalised spectrum to file
   np.savetxt("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.dat".format(pathtodat,iso,lab,mass),scaled_spectrum)

   #double check
   data = np.loadtxt("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.dat".format(pathtodat,iso,lab,mass))
   fig = plt.figure() #figsize=(8,6))
   ax = fig.add_subplot(111)
   ax.plot(scaled_spectrum["ER"],    scaled_spectrum["dRdE"])
   ax.plot(data[:,0],data[:,1])
   ax.set_xlabel("Recoil energy (keV)")
   ax.set_ylabel("dR/dE (per keV)")
   fig.savefig("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.png".format(pathtodat,iso,lab,mass))

iso = "Comb"  #131
llist = ["c{0}p=c{0}n".format(i) for i in [1,3,8]]
mlist = [m for m in [50,100,500,1000]]

for l in llist:
   for m in mlist:
      get_norm_spectra(iso, l, m)

