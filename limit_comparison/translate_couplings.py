import numpy as np
import limit_comparison.parameters as p

pi = np.pi

# Scalar-nucleon coupling from scalar-quark and scalar-gluon couplings
# eq. B.9
# Here I use cG = (c_g/Lambda)*(alpha_s/(12*pi)) to match other results more easily
def calc_cN_scalar(cu,cd,cs,cc,cb,ct,cG,data):
   mN = data["mN"]
   mu = data["mu"]
   md = data["md"]
   ms = data["ms"]
   mc = data["mc"]
   mb = data["mb"]
   mt = data["mt"]
   f_Tu = data["f_Tu"]
   f_Td = data["f_Td"]
   f_Ts = data["f_Ts"]
   f_TG = 1 - (f_Tu + f_Td + f_Ts)

   cgLinv = cG * (12*pi/data["alpha_s"])

   Xu = cu*(mN/mu)
   Xd = cd*(mN/md)
   Xs = cs*(mN/ms)
   Xc = cc*(mN/mc)
   Xb = cb*(mN/mb)
   Xt = ct*(mN/mt)
   Xg = cgLinv*mN  # i.e. cg*(mN/Lambda)
   cN = Xu*f_Tu \
      + Xd*f_Td \
      + Xs*f_Ts \
      + (2./27.)*f_TG*(Xc+Xb+Xt - Xg)
   return cN

# eq. B.16
# Here I use cG = (c_g/Lambda)*(alpha_s/(8*pi)) to match other results more easily
def calc_cN_pseudoscalar(cu,cd,cs,cc,cb,ct,cG,data):
   mN = data["mN"]
   mu = data["mu"]
   md = data["md"]
   ms = data["ms"]
   mc = data["mc"]
   mb = data["mb"]
   mt = data["mt"]
   Delu = data["Delu"]
   Deld = data["Deld"]
   Dels = data["Dels"]

   cgLinv = cG * (8*pi/data["alpha_s"])

   mhat = 1./(1./mu + 1./md + 1./ms)
   Cx = cu*(mhat/mu) \
      + cd*(mhat/md) \
      + cs*(mhat/ms) \
      + cc*(mhat/mc) \
      + cb*(mhat/mb) \
      + ct*(mhat/mt) \
      - cgLinv*mhat  # cg*(mhat/Lambda)
   Yu = (mN/mu)*(cu - Cx)*Delu
   Yd = (mN/md)*(cd - Cx)*Deld
   Ys = (mN/ms)*(cs - Cx)*Dels
   cN = Yu + Yd + Ys
   return cN 

# eq. B.18
def calc_cp_vector(cu,cd):
   return 2*cu + cd
def calc_cn_vector(cu,cd):
   return cu + 2*cd

# eq. B.20
def calc_cN_axialvector(cu,cd,cs,data):
   Delu = data["Delu"]
   Deld = data["Deld"]
   Dels = data["Dels"]

   return cu*Delu + cd*Deld + cs*Dels

# eq. B.22
def calc_cN_tensor(cu,cd,cs,data):
   delu = data["delu"]
   deld = data["deld"]
   dels = data["dels"]
   
   return cu*delu + cd*deld + cs*dels

#def calc_cN_axialtensor(cu,cd,cs,data)
# need to calculate


# Conversion factors for WIMP-nucleon coupling to WIMP-nucleon cross section
#===================================================
# i.e. if one has the coupling squared, multiply it by the following factors
# to get the cross-section

# Reduced mass
def Mred(m1,m2):
  return m1*m2/(m1+m2)

# General list of factors and functions
gsq_to_sigmaOnMr2 = []
def gsq_to_sigma(mchi,Nop):  # fold in the reduced mass
   return Mred(mchi,p.mN)**2 * gsq_to_sigmaOnMr2[Nop]

# O0 (placeholder to offset indices to match operator number)
gsq_to_sigmaOnMr2 += [None]
# O1 (SI) , g1^2 --> sigma_SI 
gsq_to_sigmaOnMr2 += [ 1 / pi * p.GeV2_to_cm2 ]  # has reduced mass factored out
# O2 (unused)
gsq_to_sigmaOnMr2 += [None]
# 03
gsq_to_sigmaOnMr2 += [None]
# O4 (SD) , g4^4 --> sigma_SD   #formula from tools paper
gsq_to_sigmaOnMr2 += [0.25 / pi * p.GeV2_to_cm2]  # has reduced mass factored out
# "Tools" paper uses '3'. Factor above is just tuned so that EFT results match Xenon225 limit. Need to check this calculation
# However, factor of 3 seems to better match the below 'sigmaSD_M6q_v2'



