import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

labels = [["c{0}p".format(i),"c{0}n".format(i),"c{0}p=c{0}n".format(i)] for i in range(1,16)]

masses = np.array([5,10,50,100,500])
mlabels = ["{0} GeV".format(m) for m in masses]


# Loop through operators
# Compute total number of events
# Compute coupling required to get 10 events

iso = "Comb"

# Loop through operators
for i,l in enumerate(labels):
   # Loop through p,n,p+n combinations
   for p,lab in enumerate(l):
      data = np.loadtxt("/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/Xe{0}/Xe{0}_{1}.dat".format(iso,lab))
      for m,mass in enumerate(masses):
         print "mWIMP={0}, {1}".format(mass,lab)
         spectrum = data[:,[0,m+1]]
         norm = 1 # should be 1 datapoint per keV
         rate0 = np.sum(spectrum[:,1])/norm
         c0 = 1
         #(rate0/rate1 = c0^2 / c1^2,  i.e. c1^2 = c0^2 * (rate1/rate0) 
         print "  Expected total events with c=1: ", rate0
         print "  Coupling required for 1   events: c={0}".format( np.sqrt(c0**2 * (1/rate0)) )
         print "  Coupling required for 10  events: c={0}".format( np.sqrt(c0**2 * (10/rate0)) )
         print "  Coupling required for 100 events: c={0}".format( np.sqrt(c0**2 * (100/rate0)) )

