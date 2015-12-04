# Conversion factor for GeV^2 to cm^2
GeV2_to_cm2 = 0.3894*10**-27

mWeak = 246.2 # GeV

# First do operators M1 and M7, easier.
# These just contribute to the SI interaction

# Constants (1307.5955 app. B + PDG)
# Quark masses just taken as PDG mean values, not entirely consistently (e.g. some MSbar (not all at same scale), some pole)
alpha_s = 0.1185 # just value at mZ for now
mN = 0.938
mu = 0.0023
md = 0.0048
ms = 0.095
mc = 1.275
mb = 4.18
mt = 173.2
fp_Tu =  0.023
fp_Td =  0.034
fN_Ts =  0.14
fn_Tu =  0.019
fn_Td =  0.041
Delp_u =  0.77
Delp_d = -0.4
DelN_s = -0.12

# Nucleon data
data_N = {  "alpha_s": alpha_s
          , "mN": mN 
          , "mu": mu
          , "md": md
          , "ms": ms
          , "mc": mc
          , "mb": mb
          , "mt": mt
          , "f_Ts": fN_Ts
          , "Dels": DelN_s
        }

# Proton data (also ignoring uncertainties and just taking one set of results)
data_p = {  "f_Tu":  fp_Tu
          , "f_Td":  fp_Td
          , "Delu":  Delp_u 
          , "Deld":  Delp_d
         }
data_p.update(data_N)

# Neutron data (as above)
data_n = {  "f_Tu":  fn_Tu
          , "f_Td":  fn_Td
          , "Delu":  Delp_d
          , "Deld":  Delp_u
         }
data_n.update(data_N)

