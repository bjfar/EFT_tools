import numpy as np

labels = [["c{0}p".format(i),"c{0}n".format(i),"c{0}p=c{0}n".format(i)] for i in range(1,16)]

masses = np.array([5,10,50,100,500])
isotopes = [128,129,130,131,132,134,136]

# original exposure, needed so we can rescale results for some different exposure 
orig_exposure = 7800. #in kilogram.days

def normalise_spectrum(in_spectrum,desired_exposure,normrate=5,normrange=(0,1000),c0=1,norm=1,verbose=True):
   # in_spectrum: (2,n) array, column 1 is recoil energy in keV with normalisation=norm
   #                           column 2 is dR/dE in (keV)^-1
   # desired_exposure: Exposure in kg.days to use for coupling computations
   # normrate: Total expected number of events within normrange desired for rescaled spectrum
   # normrange: Range of recoil energies over which to compute normrate
   # c0: value of coupling used to generated original spectrum
   # norm: e.g. 1 if datapoints are separated by 1 keV

   # Find the desired normrange endpoints in ER indices
   #startcut = np.argmin(np.abs(spectrum[:,0] - normrange[0]))
   #endcut   = np.argmax(np.abs(spectrum[:,0] - normrange[1]))
   mask = (in_spectrum[:,0] >= normrange[0]) & (in_spectrum[:,0] <= normrange[1])

   #mask = (in_spectrum[:,0] >= normrange[0]) & (in_spectrum[:,0] <= normrange[1])
   #print "sum mask: ", np.sum(mask)
   #print "sum mask: ", np.min(np.where(mask)), np.max(np.where(mask))
   #print "sum normrange: ", normrange
   #print "rate0 ({0}) = {1}".format(normrange,np.sum(in_spectrum[mask][1]))
   #print "rate0 ({0}) = {1}".format(normrange,np.sum(in_spectrum[normrange[0]:normrange[1],1]))

   if verbose: print "  Expected total events (in range {0}, with original exposure={1} kg.days) = {2}".format(normrange, orig_exposure, np.sum(in_spectrum[mask][:,1])/norm)

   exposeK = desired_exposure/orig_exposure
   
   # rescale original spectrum for the desired exposure
   spectrum = np.empty_like(in_spectrum) 
   spectrum[:] = in_spectrum #make sure copy happens...
   spectrum[:,1] = exposeK*spectrum[:,1]

   # Compute total number of events in normrange
   rate0 = np.sum(spectrum[mask][:,1])/norm
   if verbose: print "  Expected total events (in range {0}, with desired exposure ={1} kg.days) = {2}".format(normrange, desired_exposure, rate0)
   
   # Scale spectra so that there are 10 expected events in total
   #norm = 1 # should be 1 datapoint per keV

   #(rate0/rate1 = c0^2 / c1^2,  i.e. c1^2 = c0^2 * (rate1/rate0) 
   # i.e. rate1 = (c0^2 / c1^2) * rate0
   c1   = np.sqrt(c0**2 * (normrate/rate0))
   scaling1 = normrate / rate0  # = (c1**2 / c0**2) 
   #np.sum(scaled_spectrum["dRdE"][mask])
   print "test:"
   print scaling1
   print np.sum(scaling1*spectrum[mask][:,1])
   print np.sum((c1**2/c0**2)*spectrum[mask][:,1])
   #print spectrum[8:30,0]
   #print spectrum[mask,0]

   #fraction of events below 30 keV / 200 keV
   scaled_spectrum = np.array(np.zeros(spectrum.shape[0]),dtype=np.dtype([("ER",np.int16),("dRdE",np.float16)]))
   scaled_spectrum["ER"] = spectrum[:,0]
   scaled_spectrum["dRdE"] = scaling1*spectrum[:,1]
   print np.sum(scaled_spectrum["dRdE"][mask])

   #print scaled_spectrum
   inittotal= np.sum(spectrum[2:,1])
   total   = np.sum(scaled_spectrum["dRdE"][7:])
   frac30  = np.sum(scaled_spectrum["dRdE"][7:30])/total
   frac200 = np.sum(scaled_spectrum["dRdE"][7:200])/total
   if verbose: print "  Expected total events with c=1: {0}".format(inittotal)
   if verbose: print "  Expected total events with c={0}: {1}".format(c1, total)
   if verbose: print "  Expected total events with c={0} (in range {1}, with desired exposure ={2} kg.days) = {3}".format(c1, normrange, desired_exposure, np.sum(scaled_spectrum["dRdE"][mask]))
   if verbose: print "  fraction of events below 30 keV (200 keV): {0} ({1})".format(frac30,frac200)  
   if verbose: print "  Coupling (coupling^2) required for {0} events (exposure={1} kg.days): c={2}, c^2={3}".format(normrate, desired_exposure, c1  ,c1**2  )
   #plt.plot(scaled_spectrum["ER"],    scaled_spectrum["dRdE"])
   #plt.plot(scaled_spectrum["ER"], 10*scaled_spectrum["dRdE"])
   #plt.plot(scaled_spectrum["ER"],100*scaled_spectrum["dRdE"])
   #plt.xlabel("Recoil energy (keV)")
   #plt.ylabel("dR/dE (per keV)")
   #plt.show()
   return scaled_spectrum, c1

   # # output normalised spectrum to file
   # np.savetxt("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.dat".format(pathtodat,iso,lab,mass),scaled_spectrum)

   # #double check
   # data = np.loadtxt("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.dat".format(pathtodat,iso,lab,mass))
   # fig = plt.figure() #figsize=(8,6))
   # ax = fig.add_subplot(111)
   # ax.plot(scaled_spectrum["ER"],    scaled_spectrum["dRdE"])
   # ax.plot(data[:,0],data[:,1])
   # ax.set_xlabel("Recoil energy (keV)")
   # ax.set_ylabel("dR/dE (per keV)")
   # fig.savefig("{0}/Xe{1}_{2}_M={3}GeV_EinkeV_sumdRdE=1.png".format(pathtodat,iso,lab,mass))


