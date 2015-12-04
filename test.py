import numpy as np
import matplotlib.pyplot as plt

pathtodat = "/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/"

Ops = [(3,"O3","0_3"),(8,"O8","O_8")]

for Onum,Oname,Olatexname in Ops:
   datafile_basenames = ["CDMSIIGe", "CDMSIISi", "LUX", "SuperCDMS_Soudan"]
   datafiles = ["{0}_{1}.dat".format(base,Oname) for base in datafile_basenames]
   CDMSIIGe, CDMSIISi, LUX, SuperCDMS_Soudan = [np.loadtxt("{0}/{1}".format(pathtodat,fname),delimiter=",") for fname in datafiles]
      
   #Sort the curve data
   cdata = [CDMSIIGe, CDMSIISi, LUX, SuperCDMS_Soudan]
   csort = [c[c[:,0].argsort()] for c in cdata]

   print csort[0].T

   plt.plot(*csort[0].T)
   plt.show()
   quit()
