import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

#userel = True # Switch from non-relativistic to relativistic operators
userel = False

if userel:
   print "Assuming input is data from RELATIVISTIC operators"
   nops = 20
   tag = "R"
else:
   print "Assuming input is data from NON-relativistic operators"
   nops = 15
   tag = "NR"

labels = [["c{0}p".format(i),"c{0}n".format(i),"c{0}p=c{0}n".format(i)] for i in range(1,nops+1)]

#masses = np.array([5,6,10,30,50,100,300,500,700,800,1000,2000,3000,5000]) #All currently available masses
masses = np.array([6,10,50,100,500,1000,5000])
mlabels = ["{0} GeV".format(m) for m in masses]
colors = ["blue","red","limegreen","m","c","gold","gray","blue","red","limegreen","m","c","gold","gray","blue","red","limegreen","m","c","gold","gray"] 

# Want to compare values against the maximum kinetic energy possessed by any WIMP in the halo:

c = 2.998*10**8 #m/s
c2 = c**2
mNucleon = 0.987 #GeV/c^2
v0 = 220 * 1000 #m/s (mean WIMP speed in galactic rest frame)
vesc = 550 * 1000 #m/s (halo escape velocity)
mT = 132*mNucleon;
muT = masses*mT/(masses + mT)
ERmaxesc = (2 * (muT/c2)**2 * vesc**2) / (mT/c2) * 10**6
ERmax0   = (2 * (muT/c2)**2 * v0**2) / (mT/c2) * 10**6

#Isotope %      spin
#124Xe   0.09 	# 0 	
#126Xe   0.09 	# 0 	
#128Xe   1.92 	# 0 	
#129Xe   26.44 	# 1/2 
#130Xe   4.08   # 0   
#131Xe   21.18	# 3/2 
#132Xe   26.89	# 0 	
#134Xe   10.44 	# 0 	
#136Xe   8.87  	# 0

XeAbun = {
  124: 0.0009, 
  126: 0.0009, 
  128: 0.0192, 
  129: 0.2620,  #ran
  130: 0.0408, 
  131: 0.2180,  #ran
  132: 0.2689,
  134: 0.1044,
  136: 0.0887,
}

# Define subset of isotopes to examine and reweight abundances to
# compensate for missing isotopes.
isotopes = [128,129,130,131,132,134,136]
regencomb = True # Regenerate the combined data?
dosep = True # Generate plots for individual isotopes 

tot = np.sum([XeAbun[iso] for iso in isotopes])
relw = np.array([XeAbun[iso]/tot for iso in isotopes])

def makeplots(axes,data,labels):
   for p,lab in enumerate(labels):
      ax = axes[p]
      print "makeplots:", p, lab, data.shape
      for m,mass in enumerate(masses):
         ymax = np.max(data[:,:,m+1])
         ax.plot(data[p,:,0],data[p,:,m+1]/ymax,label=mass,color=colors[m])
         ax.axvline(ERmax0[m]  ,color=colors[m],linestyle="dashed")
         ax.axvline(ERmaxesc[m],color=colors[m],linestyle="dotted")
      ax.set_xlim(0,600)
      ax.set_title(lab,fontsize=8, y=0.8)
      leg = ax.legend(loc=1, frameon=False, framealpha=0,prop={'size':6})

if regencomb:
   print "Regenerating combined recoil spectra"
   # Loop through operators
   for i,l in enumerate(labels):
      # Loop through p,n,p+n combinations
      for p,lab in enumerate(l):
         dataiso = []
         # Loop through isotopes
         for n,(iso,w) in enumerate(zip(isotopes,relw)):
            print iso, w, lab
            datai = np.loadtxt("/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/Xe{0}/Xe{0}_{1}_{2}.dat".format(iso,tag,lab))
            print "Data array shape (should stay consistent): {0}".format(datai.shape)
            dataiso += [datai]
         # dataiso should now have shape (n_isotopes, n_ER, n_DMmasses+1)
         # Combine spectra weighted by isotope abundance
         dataiso = np.array(dataiso) # should have shape (n_isotopes, n_DMmasses+1, n_ER)
         datacomb = np.sum(dataiso[:,1:]*relw[:,np.newaxis,np.newaxis],axis=0) #sum over isotopes (dim 0)
         np.savetxt("/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/XeComb/XeComb_{0}_{1}.dat".format(tag,lab),datacomb)


print "Beginning plot generation..."

# Loop through isotopes (plus combined)
todolist = ["Comb"]
if dosep: todolist += isotopes
for iso in todolist:
   fig, axes = plt.subplots(nops,3,figsize=(12,2*nops),sharey=True)
   # Loop through operators
   for i,l in enumerate(labels):
      data = []
      # Loop through p,n,p+n combinations
      for p,lab in enumerate(l):
         datap = np.loadtxt("/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/Xe{0}/Xe{0}_{1}_{2}.dat".format(iso,tag,lab))
         # datap has shape (n_ER, n_DMmasses+1)
         data += [datap]
   
      #"data" should have shape (n_labels (per operator), n_ER, n_DMmasses+1)
   
      data = np.array(data)
      data[np.where(data<0)]=0 #remove negative values
      print "data shape: ", data.shape
      makeplots(axes[i,:],data,labels[i])

   fig.subplots_adjust(wspace=0.01)
   name = "Xe{0}_{1}_EFTcoeffrecoilspectra.png".format(iso,tag)
   print "Saving plot {0}...".format(name)
   fig.savefig(name,bbox_inches='tight')

