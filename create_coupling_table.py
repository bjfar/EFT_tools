import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import normalise_spectra as ns
from scipy.interpolate import interp1d

pathtodat = "EFTcoeffplotdata/"

#labels = [["c{0}p".format(i),"c{0}n".format(i),"c{0}p=c{0}n".format(i)] for i in range(1,16)]
operators = ["c{0}p=c{0}n".format(i) for i in range(1,16)]

allmasses = np.array([5,6,10,30,50,100,300,500,700,800,1000,2000,3000,5000]) 
#masses = np.array([6,10,50,100,300,1000])
masses = allmasses
colors = ["blue","k","darkred","limegreen","m","darkblue","goldenrod","gray","orange","darkgreen","c","green","k","purple","red"]
isotopes = [128,129,130,131,132,134,136]

# Exposure to use for normalisation
exposure = 225*34 * 0.5  # 50% exposure to roughly account for cuts/acceptance
#exposure = 85*118   # LUX

iso = "Comb"

# Assumed number of events observed
Nevents = 5

# Basic wiki table structure:
# {|
# |-
# | A || B
# |-
# | C || D
# |}
# Will generate table of this format.
def format_table(table,colwidth,extracols):
    newtable = ""

    firstrow = "|-\n|     | "
    for j,m in enumerate(masses):  
       firstrow += "|| {0:{1}d} GeV ".format(m,colwidth)
       for ecol in extracols:
          if(ecol[0]==m):
             firstrow += "|| {0:{1}d} GeV (CDSM) ".format(m,colwidth)
    newtable += firstrow + "|  |\n"

    for i,line in enumerate(table):
        row = ""
        if np.mod(i,3)==0:
           col1 = "|-\n| O{0:<{1}} | 8-30 keV  ".format(i/3+1,2)
           lastcol = "| SuperCDMS |"
        elif np.mod(i,3)==1:
           col1 = "|  {0:<{1}} | 8-240 keV  ".format("",2)
           lastcol = "| CDMSIIGe  |"
        else:
           col1 = "|  {0:<{1}} | 8-1000 keV ".format("",2)
           lastcol = "| CDMSIISi  |"
        row += col1
        for j,m in enumerate(masses):
           if line[j]<0:
              row += "|| **{0:{1}.2e}**".format(-line[j],colwidth)
           else:
              row += "|| {0:{1}.2e}".format(line[j],colwidth)
           for ecol in extracols:
              if(ecol[0]==m):
                 print ecol, i
                 row += "|| {0:{1}.2e}".format(ecol[1][np.mod(i,3)][i//3],colwidth)           
        #row = col1+" | " + " | ".join("{0:{1}.2e}".format(x,colwidth) for x in line) + " |"           
        newtable += row + lastcol + "\n"
       #if np.mod(i,3)==2:
           #newtable += ["____________________________________________________________________________________________"]
        

    #newtable += ["{0:{1}}".format("",3)+" | " + " | ".join("{0:{1}d}".format(x,9) for x in masses) + " |"]
    return newtable

table = []
for O in operators:
   rowa = []
   rowb = []
   rowc = []
   for m in masses:
      print "Extracting: ",iso, O, "mWIMP=",m
      data = np.loadtxt("{0}/Xe{1}/Xe{1}_{2}.dat".format(pathtodat,iso,O))
      print "Normalising spectrum for Xe iso={0}, operator={1}, mWIMP={2}".format(iso,O,m)
      index = np.where(allmasses==m)[0][0]+1 #first column is ER
      print index
      spectrum = data[:,[0,index]] 
      scaled_spectrum, c1a = ns.normalise_spectrum(spectrum,exposure,normrate=Nevents,normrange=(8,30),verbose=True)
      scaled_spectrum, c1b = ns.normalise_spectrum(spectrum,exposure,normrate=Nevents,normrange=(8,240),verbose=True)
      scaled_spectrum, c1c = ns.normalise_spectrum(spectrum,exposure,normrate=Nevents,normrange=(8,1000),verbose=True)
      c1a2 = c1a**2
      c1b2 = c1b**2
      c1c2 = c1c**2

      # Sanity check:
      # should have c1a2 > c1b2 > c1c2
      if (c1c2 > c1b2) or (c1c2 > c1a2):
         raise ValueError("Relative coupling sizes make no sense! Should have c1a^2 > c1b^2 > c1c^2, but have instead {0},{1},{2}".format(c1a2,c1b2,c1c2))

      if c1a**2>2*c1b**2: c1b2 = -c1b**2
      if c1b**2>1.2*c1c**2: c1c2 = -c1c**2
      rowa+=[c1a2]
      rowb+=[c1b2]
      rowc+=[c1c2]
   table += [rowa]
   table += [rowb]
   table += [rowc]

table = np.array(table)

# Limits from CDMS paper (arxiv:)
lims = {}
lims["SuperCDMS_10GeV"] =[8.98e-5,np.nan,3.14e4,8.77e1,6.34e5,4.54e8 ,8.44e7,4.30e2,1.95e5,9.22e4,5.13e-1,1.03e2,4.28e8 ,5.00e11,1.32e8]
lims["SuperCDMS_300GeV"]=[np.nan for i in range(15)]

lims["CDMSIIGe_10GeV"]  =[2.00e-3,np.nan,2.24e5,2.05e3,9.18e6,3.30e9,2.51e9,1.16e4,2.48e6,1.11e6,6.15e0,1.21e3,3.06e9,8.20e12,5.65e8]
lims["CDMSIIGe_300GeV"] =[8.42e-6,np.nan,2.66e1,1.10e1,4.04e3,4.50e5,1.12e7,2.67e1,3.87e3,9.08e2,5.46e-3,8.70e-1,3.56e5,8.46e9,1.10e4]
# What even are these? Numbers seem to be totally wrong...
#lims["CDMSIIGe_10GeV"]  =[1.19e-3,np.nan,1.06e5,1.24e3,5.30e6,1.55e9 ,1.76e9,7.68e3,1.32e6,5.83e5,3.23e0 ,6.33e2,1.44e9 ,4.91e12,2.76e8]
#lims["CDMSIIGe_300GeV"] =[1.13e-5,np.nan,3.08e1,1.53e1,4.82e3,5.21e5 ,1.62e7,3.51e1,4.84e3,1.09e3,6.59e-3,1.04e0,4.12e5 ,1.06e10,1.26e4]

lims["CDMSIISi_10GeV"]  =[3.06e-3,np.nan,8.59e5,3.94e3,2.67e7,2.44e10,3.19e9,1.70e4,9.17e6,4.34e6,1.86e1 ,2.45e3,2.50e13,2.64e13,4.44e9]
lims["CDMSIISi_300GeV"] =[7.73e-4,np.nan,1.37e4,1.02e3,1.55e6,3.70e8 ,9.29e8,3.49e3,7.34e5,2.86e5,1.34e0 ,1.68e2,1.36e12,1.72e12,1.48e7]
extracols=[(10,(lims["SuperCDMS_10GeV"],lims["CDMSIIGe_10GeV"],lims["CDMSIISi_10GeV"])),
           (300,(lims["SuperCDMS_300GeV"],lims["CDMSIIGe_300GeV"],lims["CDMSIISi_300GeV"]))]

newtable = format_table(table,9,extracols)
#newtable += ["{0:{1}}".format("",3)+" | " + " | ".join("{0:{1}d}".format(x,9) for x in masses) + " |"]

fname = "{0}/coupling_table.txt".format(pathtodat)
print "Writing coupling table to {0}".format(fname)
with open(fname,"w") as f:
   f.write(newtable)

# All operators on one plot
plotOs = [(1,3,4,5),(6,7,8,9,10),(11,12,13,14,15)]
#for i,row in enumerate(table):

for k,plot in enumerate(plotOs):
   fig = plt.figure(figsize=(6,4))
   ax = fig.add_subplot(111)
   for Onum in plot:

      finemasses = np.exp(np.arange(np.log(masses[0]),np.log(masses[-1]),np.log(1.01)))
      kind = "cubic"
      f1 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3])),   kind=kind, bounds_error=False)
      f2 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3+1])), kind=kind, bounds_error=False)
      f3 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3+2])), kind=kind, bounds_error=False)
 
      ax.set_xlabel(r"$\mathrm{WIMP\,mass\,[GeV]}$")
      ax.set_ylabel(r"$c^2 * (m^4_\mathrm{weak})$")
      ax.set_xscale("log")
      ax.set_yscale("log")

      ax.set_xlim(3,1000)
      #ax.set_ylim(1e-3,1e9)
      #ax.grid(True)

      #row = np.abs(table[(Onum-1)*3])
      #ax.plot(masses, row, c=colors[Onum-1])

      ax.plot(finemasses, np.exp(f1(np.log(finemasses))), label="$O_{{{0}}}$".format(Onum),  c=colors[Onum-1], lw=1.5 )
      ax.plot(finemasses, np.exp(f2(np.log(finemasses))), c=colors[Onum-1], ls="dashed")
      ax.plot(finemasses, np.exp(f3(np.log(finemasses))), c=colors[Onum-1], ls="dotted")
       
      #ax.text(0.7*masses[0],table[(Onum-1)*3][0],"O{0}".format(Onum))

      leg = ax.legend(frameon=False, framealpha=0, prop={'size':10}, ncol=2)
      for legobj in leg.legendHandles:
        legobj.set_linewidth(2.0)
      plt.tight_layout()

      fig.savefig("{0}/N={1}_ops_group{2}.png".format(pathtodat,Nevents,k))
      print "created {0}/N={1}_ops_group{2}.png".format(pathtodat,Nevents,k)

# Nicer plots

# Compare to CDMS:
Ops = [(3,"O3","O_3"),(8,"O8","O_8")]
#Ops = [(8,"O8","O_8")]

for Onum,Oname,Olatexname in Ops:
   datafile_basenames = ["CDMSIIGe", "CDMSIISi", "LUX", "SuperCDMS_Soudan"]
   datafiles = ["{0}_{1}.dat".format(base,Oname) for base in datafile_basenames]
   curvedata = [np.loadtxt("{0}/{1}".format(pathtodat,fname)) for fname in datafiles]
      
   # Sort the curve data in increasing x value
   CDMSIIGe, CDMSIISi, LUX, SuperCDMS_Soudan = [c[c[:,0].argsort()] for c in curvedata]

   finemasses = np.exp(np.arange(np.log(masses[0]),np.log(masses[-1]),np.log(1.01)))
   kind = "cubic"
   f1 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3])),   kind=kind, bounds_error=False)
   f2 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3+1])), kind=kind, bounds_error=False)
   f3 = interp1d(np.log(masses), np.log(np.abs(table[(Onum-1)*3+2])), kind=kind, bounds_error=False)
 
   fig = plt.figure(figsize=(6,4))
   ax = fig.add_subplot(111)
   ax.plot(*CDMSIIGe.T, label="CDMSII Ge 90% CL", c="g")
   ax.plot(*CDMSIISi.T, label="CDMSII Si 90% CL", c="grey")
   ax.plot(*SuperCDMS_Soudan.T, label="SuperCDMS Soudan 90% CL", c="m")
   ax.plot(*LUX.T,      label="LUX 90% CL (est. by CDMS)", c="b")
   ax.plot(finemasses, np.exp(f1(np.log(finemasses))), label="Xenon100 N={0} in 8-30 keV".format(Nevents),  c="r", lw=2 )
   ax.plot(finemasses, np.exp(f2(np.log(finemasses))), label="Xenon100 N={0} in 8-240 keV".format(Nevents), c="r", ls="dashed")
   ax.plot(finemasses, np.exp(f3(np.log(finemasses))), label="Xenon100 N={0} in 8-1000 keV".format(Nevents), c="r", ls="dotted")
   ax.set_xlabel(r"$\mathrm{WIMP\,mass\,[GeV]}$")
   ax.set_ylabel(r"$c^2 * (m^4_\mathrm{weak})$")
   ax.set_xscale("log")
   ax.set_yscale("log")
   ax.set_xlim(3,1000)
   ax.set_ylim(1e-1,1e8)
   #ax.grid(True)
   ax.set_title("${0}$".format(Olatexname))
   ax.legend(frameon=False, framealpha=0, prop={'size':10})
   plt.tight_layout()
   fig.savefig("{0}/N={1}_{2}.png".format(pathtodat,Nevents,Oname))
   print "created {0}/N={1}_{2}.png".format(pathtodat,Nevents,Oname)
#plt.show()

#newtable += []
   
