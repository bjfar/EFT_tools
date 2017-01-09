#!/usr/local/bin/MathKernel -script

(* Compute recoil spectra for a specified isotope of Xenon *)

(*Args:*)
(*1 - Absolute path to "v6/dmformfactor-V6.m"
      or to "v4/dmformfactor-V4-IDM.m" for inelastic mode (package from Chang et. al. 1409.0536) *)
(*2 - Relative output path*)
(*3 - Xenon isotope (e.g. "131")*)
(*4 - flag to use relativistic operator coefficients*)
(*5 - flag to switch form factor treatment to traditional Helm form factors rather than derived from nuclear density matrices *) 
(*6 - flag to switch to inelastic scattering mode *)

(* Make output more readable in the terminal *)
(* SetOptions["stdout", FormatType->InputForm] *)

(* Turn on handler to abort evaluation if anything weird happens *)
messageHandler = If[Last[#], Print["An error has occurred in 'recoil_generator.m'; aborting evaluation..."] Exit[1];]&
Internal`AddHandler["Message", messageHandler]
(* ..but need to suppress expected warnings *)
Off[Part::partd]
Off[Simplify::time]

args = Join[$CommandLine[[4;;]],Table[{},{i,1,10}]];
(* Print["Running with command line arguments: "<>ToString[args]]; *)
outpath = Directory[]<>"/"<>args[[2]];
Print["Setting output directory to: "<>outpath];
If[!DirectoryQ[outpath],
  Print["\nOutput directory "<>outpath<>" does not exist! Please create it and try again.\n"];
  Exit[1];
];

path = args[[1]];
isotope = ToExpression[args[[3]]];
Otype = args[[4]]
checkrel = SameQ[args[[5]],"check"]; (* Compute coefficients of relativistic operators using their non-relativistic reductions (to check that I understand the mapping correctly *)
helm = SameQ[args[[6]],"UseHelm"]; (* Use "traditional" Helm form factors rather than ones derived from nuclear density matrices *)
idm = SameQ[args[[7]],"IDM"]; (* Alter scattering kinematics to inelastic; code from 1409.0536; note, must also alter package path*)

(* For Notebook testing, delete everything above this and uncomment the following: *)
(*
outpath = "recoil_spectrum_tables";
path = "/home/farmer/mathematica/DMFormFactor_13086288";
isotope = 131;
Otype = "NR";
checkrel = False; (* Compute coefficients of relativistic operators using their non-relativistic reductions (to check that I understand the mapping correctly *)
helm = False; (* Use "traditional" Helm form factors rather than ones derived from nuclear density matrices *)
*)

userel = SameQ[Otype,"R"]; (* Switch to using coefficients of relativistic operators (spin 1/2 WIMP only) *)
Print["Looking for DMFormFactor in directory: "<>path]
If[!DirectoryQ[path],
  Print["\nSupplied DMFormFactor directory "<>path<>" does not exist! Please check it for typos and try again\n"];
  Exit[1];
];
Check[
  If[idm,
    Import[path<>"/v4/dmformfactor-V4-IDM.m"];
    ,
    Import[path<>"/v6/dmformfactor-V6.m"]; 
  ]
 ,
  If[idm,
    Print["...Error loading DMFormFactor package. Please check the path (and perhaps version; you have selected INELASTIC mode so this script requires the modified v4 from 1409.0536) and try again"];
    ,
    Print["...Error loading DMFormFactor package. Please check the path (and perhaps version; this script requires v6) and try again"];
  ]
  Exit[1];
]
Print["...Found."];

(* Check that chosen Xenon isotope is valid *)
validiso = {128,129,130,131,132,134,136};
If[Count[validiso,isotope]==0, 
  Print["\nRequested Xenon isotope '"<>ToString[isotope]<>"' is invalid! Aborting...\n"]; 
  Exit[1];
];

(* Helper functions *)
SetAttributes[DoSilent, HoldAll];
DoSilent[expr_] := Block[{Print = Null &}, expr];

(* WIMP masses for which we want recoil spectra *)
(*m\[Chi]={5,10,50}; (* GeV *)*)

(* This one used a lot up until R vs NR bug found, now switching to limited list below*)
(*m\[Chi]={3,4,5,6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,100,200,300,500,700,800,1000,2000}; (* GeV *)*)
m\[Chi]={5,6,7,8,9,10,16,20,30,40,50,70,100,300,500,800,1000,3000,5000,8000,10000};

(*m\[Chi]={3,5,6,7,8,9,10,15,20,30,50,100,300,500,700,800,1000,2000,3000,5000}; (* GeV *) *)
(*m\[Chi]={5,50,500,5000}; (* GeV *) *)
spin=1/2 (*leave hardcoded for now*)

(* Specifiy recoil energies at which to evaluate dR/dE, in keV *)
(* Two versions: one for low PE signal region, one for high PE
   Goes well beyond signal region boundaries to account for detector response (smearing) *)
Earrlow =Table[N[ER*10^-6],{ER,0.05,100,0.05}];
Earrhigh=Table[N[ER*10^-6],{ER,0.5,1000,0.5}];

(* For INELASTIC scattering; select delm (mass splitting) values to use (in keV)*)
delmlist={0,5,10,30,70,100,150,200,250,300,400,500,700,1000};

(* Exposure to use for rate calculations. Can be easily scaled to something different in later analysis steps. *) 
benchmarkExposure = 7800 KilogramDay;

Print["**************************************"    ]
Print[" Generating recoil spectra for Xe"<>ToString[isotope]<>"..."]
Print[" ...for spin "<>ToString[spin,InputForm]<>" WIMP"    ]
Print[" ...for masses (GeV) "<>ToString[m\[Chi]]  ]
Print[" ...for operator type "<>Otype             ]
If[helm,
  Print[" ...with traditional Helm form factors"    ]]
If[idm,
  Print[" ...using INELASTIC WIMP scattering kinematics"]
  Print[" ...for mass splittings (keV) "<>ToString[delmlist]  ]
]
Print["**************************************"    ]

(* Model Setup *)
SetJChi[1/2](* WIMP spin *)
SetMChi[MWIMP GeV] (* WIMP mass in GeV *)
densmatrix="default";  (* nuclear density matrix (using default density matrices) *)
SetHelm[helm]; (* If true use "traditional" Helm form factors rather than ones derived from nuclear density matrices *)
bFM="default"; (* oscillation parameter (default: approximate formulae used) *)
IsotopeA=isotope;
Print["Setting A="<>ToString[IsotopeA]]
SetIsotope[54,IsotopeA,bFM,densmatrix] (*[Z,A,bFM,filename]*)
ZeroCoeffs[];(*Reset all coefficients*)
If[userel && !checkrel,
   (* Relativistic operator coefficients *)
   Print["Using RELATIVISTIC operator coefficients"];
   DoSilent[
     Do[
       SetCoeffsRel[i,cp[[i]],"p"];
       SetCoeffsRel[i,cn[[i]],"n"];
       ,{i,1,20}];
   ];
   Nops=20;
   Tag="R";
   , (*else*)
   (* Non-relativistic operator coefficients *)
   If[userel && checkrel,
     Print["Checking RELATIVISITIC operator coefficients using their non-relativistic reductions"];
     Tag="R-check";
     Nops=20;
     ,
     Print["Using NON-relativistic operator coefficients"];
     Tag="NR";
     Nops=15;
   ]
   DoSilent[
     Do[
       SetCoeffsNonrel[i,cp[[i]],"p"];
       SetCoeffsNonrel[i,cn[[i]],"n"];
       ,{i,1,15}];
   ];
]

If[helm,
  Print["WARNING! Using traditional Helm form factors instead of the more advanced nuclear density matrix treatment!"];
  Tag = Tag <> "_HelmFF";
]

(* Parameters *)
mNucleon=0.938 GeV;
NT=1/(IsotopeA* mNucleon);(*Number of target nuclei per detector mass*)
Centimeter=(10^13 Femtometer);
rhoDM=0.3GeV/Centimeter^3;(*Local dark matter density*)
(*ve=232 KilometerPerSecond;(*Earth's velocity in galactic rest frame*)*)
ve=220 KilometerPerSecond; (* from XEPHYR *)
v0=220 KilometerPerSecond;(*Mean WIMP speed in galactic rest frame*)
vesc=544 KilometerPerSecond;
SetHALO["MBcutoff"];

(* Event rate in symbolic form *)
If[idm,
   Print["Computing symbolic event rate for INELASTIC scattering..."];
   ,
   Print["Computing symbolic event rate for elastic scattering..."];
]
DoSilent[
  If[idm,
      dRdE=EventRateInelastic[NT,rhoDM,delm 10^-6 GeV,qGeV,ve,v0,vesc];
      ,
      dRdE=EventRate[NT,rhoDM,qGeV,ve,v0,vesc]; (* as function of q, not ER *)
    ]
];


(* Produce recoil spectra for several WIMP masses *)
Print["Computing auxilliary parameters..."];
mT = IsotopeA*mNucleon;
\[Mu]T=Table[(m\[Chi][[i]] GeV mT)/(m\[Chi][[i]] GeV+mT),{i,1,Length[m\[Chi]]}];
ERmax =Table[ 2(((\[Mu]T[[i]]^2) (v^2) )/mT)/GeV/.v->ve,{i,1,Length[m\[Chi]]}];(*just using earth speed to get rough number, and get rid of GeV units for plotting*)
ERmaxesc = Table[2(((\[Mu]T[[i]]^2) (v^2) )/mT)/GeV/.v->vesc,{i,1,Length[m\[Chi]]}];
(*generalfunc= Table[2500KilogramDay*fEdRdE[ER,m\[Chi]\[LeftDoubleBracket]i\[RightDoubleBracket],mT/GeV],{i,1,Length[m\[Chi]]}];*)


(* General function to compute recoil spectra (couplings not yet replaced with values) *)
ruleQtoE1 = qGeV->Sqrt[2*mTv*ER];
ruleQtoE2 = mTv->mT/GeV;
Print["Transforming general analytic recoil spectra function..."]
generalfunc= benchmarkExposure*dRdE GeV /. ruleQtoE1 /. ruleQtoE2;

curvenames=Table["m\[Chi]="<>ToString[m\[Chi][[i]]],{i,1,Length[m\[Chi]]}];
If[idm,
  idmnames=Table["delm="<>ToString[delmlist[[i]]],{i,1,Length[delmlist]}];
]

(* Replacement rules for setting the couplings one by one *)
coefflist=Table[{cp[[i]],cn[[i]]},{i,1,20}]//Flatten;
coefflistnoiso=Table[c[[i]],{i,1,20}];
resttozero=Table[{cp[[i]]->0,cn[[i]]->0,c[[i]]->0},{i,1,20}]//Flatten;

(* Replacement rules for recreating relativistic operators by setting their non-relativistic reductions correctly (see table 1 of documentation-standalone.pdf) *)
(* Result should be equivalent to setting the various d's to 1 *)

(* Note: "The coefficients d_j are dimensionless, by inserting appropriate powers of the user-defined scale m_M (which usually would be known from the context of the theory), set by the user function SetMM. This scale is set by default to be m M = m V ≡ 246.2 GeV" *)
(* So following the above, I will use the default mM=246.2 for now *)
mM = 246.2 GeV;
mX = MWIMP GeV;
mN = mNucleon;
reductionrules=
{
  (*1*) {c[1]->1},
  (*2*) {c[10]->1},
  (*3*) {c[11]->-mN/mX},
  (*4*) {c[6]-> -mN/mX},
  (*5*) {c[1]->1},
  (*6*) {c[1]-> qGeV^2/(2*mN*mM), 
         c[3]->-2*mN/mM, 
         c[4]-> 2*mN^2/(mM*mX) * qGeV^2/mN^2,  
         c[6]->-2*mN^2/(mM*mX)},
  (*7*) {c[7]->-2,
         c[9]->2*mN/mX},
  (*8*) {c[10]->2*mN/mM},
  (*9*) {c[1]->-qGeV^2/(2*mX*mM),
         c[5]->2*mN/mM,
         c[4]->-2*mN/mM * qGeV^2/mN^2,
         c[6]-> 2*mN/mM},
  (*10*) {c[4]->4*qGeV^2/mM^2,
          c[6]->-mN^2/mM^2},
  (*11*) {c[9]->4*mN/mM},
  (*12*) {c[10]->-mN/mX * qGeV^2/mM^2,
          c[12]->-4*qGeV^2/mM^2,
          c[15]->-4*mN^2/mM^2},
  (*13*) {c[8]->2,
          c[9]->2},
  (*14*) {c[9]->-4*mN/mM},
  (*15*) {c[4]->-4},
  (*16*) {c[13]->4*mN/mM},
  (*17*) {c[11]->2*mN/mM},
  (*18*) {c[11]->qGeV^2/mM^2,
          c[15]->4*mN^2/mM^2},
  (*19*) {c[14]->-4*mN/mM},
  (*20*) {c[6]->4*mN^2/mM^2}
};

(*Get rid of unsimplified GeV units*)
Print["Simplifying out GeV units..."]
noGeVfunc = generalfunc // Simplify[#, { MWIMP \[Element] Reals, MWIMP > 0, delm \[Element] Reals, delm >= 0, GeV \[Element] Reals, GeV>0}, TimeConstraint -> 0.01] &;

Print["Creating generator tables for numerical results..."]
(* Create tables of coupling values to use *)
zeros    = Table[0,{i,1,Nops}];
zerosAll = Table[zeros,{i,1,Nops}];
cpsAll   = Table[tmp = zeros; tmp[[i]] = 1; tmp, {i, 1, Nops}];
cnsAll   = Table[tmp = zeros; tmp[[i]] = 1; tmp, {i, 1, Nops}];

(* Three sets of coupling combinations *)
(*
cCombs = {{cpsAll,zerosAll},{zerosAll,cnsAll},{cpsAll,cnsAll}};
Temporarily turning off all but isoscalar case*)
cCombs = {{cpsAll,cnsAll}};
NcCombs = 1;

(* plot titles *)
(*titles=Table[{"c"<>ToString[i]<>"p","c"<>ToString[i]<>"n","c"<>ToString[i]<>"p=c"<>ToString[i]<>"n"},{i,1,Nops}];*)
titles=Table[{"c"<>ToString[i]<>"p=c"<>ToString[i]<>"n"},{i,1,Nops}]

(* Apply the analytic expressions to the chosen recoil energies and obtain final numerical tables *)
getnum[EarrIN_] := Module[{Earr=EarrIN,out,i,m,p,d,cs,tmp,tmp2,tmpf,tmpf2,tmpf3,tmpf4},
  Print["Computing numerical results at "<>ToString[Length[Earr]]<>" recoil energies between "<>ToString[FortranForm[First[Earr]]]<>" and "<>ToString[FortranForm[Last[Earr]]]<>" GeV"];
  out=
   Table[
     WriteString["stdout", "  Evaluating operator="<>ToString[i]<>"                                ", "\n"];
     Table[
        WriteString["stdout", " Applying coefficient replacements...                     ", "\r"];
        tmpf = noGeVfunc /. {cp -> cs[[1]][[i]], cn -> cs[[2]][[i]]};
        WriteString["stdout", " Checking for equality with zero...                      ", "\r"];
        tmp = If[AllTrue[tmpf/.{delm -> 1, MWIMP -> 50, ER -> Earr},PossibleZeroQ], (*happens for some operators for some isotopes*)
          tmp2 = Table[0.*Earr, {m, m\[Chi]}]; (* set dRdE to zero *)
          If[idm,
             Table[Transpose[Join[{SetPrecision[Earr*10^6,8]},tmp2]],{d,delmlist}]
             ,
             Transpose[Join[{SetPrecision[Earr*10^6,8]},tmp2]]
          ]
          ,
          If[idm,
            Table[
              WriteString["stdout", " Applying delm = "<>ToString[d]<>" replacement...                 ", "\r"];
              tmpf2 = 10^-6 * tmpf /. {delm -> d};
              tmp2 = Table[
                WriteString["stdout", " Applying WIMP mass = "<>ToString[m]<>" replacement...                     ", "\r"];
                tmpf3 = tmpf2 /. {MWIMP -> m};
                WriteString["stdout", " Applying ER replacements...                      ", "\r"];
                tmpf4 = tmpf3 /. {ER -> Earr}
                , {m, m\[Chi]}];
              WriteString["stdout", " Joining with ER values...                  ", "\r"];
              Transpose[Join[{SetPrecision[Earr*10^6,8]},tmp2]]
            ,{d, delmlist}]
            ,
            tmp2 = Table[10^-6 * tmpf /. {MWIMP -> m, ER -> Earr}, {m, m\[Chi]}];
            Transpose[Join[{SetPrecision[Earr*10^6,8]},tmp2]]
          ]
        ]
     ,{cs,cCombs}]
   ,{i,1,Nops}];
  Print["  Complete!                                "];
  out (* "return" value *)
]

(* Function to convert numbers to Fortran-like format
from: http://mathematica.stackexchange.com/a/19547/32644 *)
f77Eform[x_?NumericQ, fw_Integer, ndig_Integer] := Module[{sig, s, p, ps},
        {s, p} = MantissaExponent[x];
     {sig, ps} = {ToString[Round[10^ndig Abs[s]]], ToString[Abs[p]]};
      StringJoin @@ Join[
                   Table[" ", {fw - ndig - 7}],
                   If[x < 0, "-", " "], {"0."}, {sig},
                   Table["0", {ndig - StringLength[sig]}], {"E"}, 
                   If[p < 0, {"-"}, {"+"}],
                   Table["0", {2 - StringLength[ps]}], {ps}]]


(* Faster write to disk than Export
from: http://mathematica.stackexchange.com/a/35375/32644 *)
writeYourCSV[file_String, list_List?MatrixQ] := 
 With[{str = OpenWrite[file, PageWidth -> Infinity], 
   len = Length[list[[1]]]}, 
  Scan[Write[str, 
     Sequence @@ (Flatten[
         Table[{FortranForm[#[[i]]], OutputForm["\t"]}, {i, 
           len - 1}]])~Join~{FortranForm[#[[len]]]}] &, list]; 
  Close[str];]

(* Save data tables to file *)
savedata[dataIN_,SRtagIN_] := Module[{data=dataIN,SRtag=SRtagIN,i,p},
  extendedoutpath=outpath<>"/Xe"<>ToString[IsotopeA];
  If[!DirectoryQ[extendedoutpath], CreateDirectory[extendedoutpath]]
  Do[
    Do[
      If[idm,
        (* old way
        Do[
           filename=extendedoutpath<>"/IDM_Xe"<>ToString[IsotopeA]<>"_"<>Tag<>"_"<>SRtag<>"_"<>titles[[i]][[p]]<>"_delm"<>ToString[delmlist[[d]]]<>".dat";
           Print["  Writing file: "<>filename];
           (*Export[filename, SetPrecision[data[[i]][[p]],8], "Table"];*)
           writeYourCSV[filename,data[[i]][[p]][[d]]];
        ,{d,1,Length[delmlist]}]
        *)
        (* Now using HDF5 format *)
        filename=extendedoutpath<>"/IDM_Xe"<>ToString[IsotopeA]<>"_"<>Tag<>"_"<>SRtag<>"_"<>titles[[i]][[p]]<>".h5";
        Print["  Writing file: "<>filename];
        Export[filename,data[[i]][[p]],{"Datasets", Table["delm"<>ToString[d],{d,delmlist}]}] 
        ,
        filename=extendedoutpath<>"/Xe"<>ToString[IsotopeA]<>"_"<>Tag<>"_"<>SRtag<>"_"<>titles[[i]][[p]]<>".h5";
        Print["  Writing file: "<>filename];
        (*Export[filename, SetPrecision[data[[i]][[p]],8], "Table"];*)
        (*writeYourCSV[filename,data[[i]][[p]]];*)
        Export[filename,data[[i]][[p]]] 
       ]
    ,{i,1,Nops}]
  ,{p,1,NcCombs}];
]
plotdatalow  = getnum[Earrlow];
savedata[plotdatalow,"lowE"];
plotdatahigh = getnum[Earrhigh];
savedata[plotdatahigh,"highE"];

(* Fix this up later
(* Report on which operators vanish for this isotope *)
Do[
  Do[
    Print[{"Vanishes? ",0==dRdE GeV//.rules[[i]][[p]]/.resttozero/.MWIMP->500/.ruleQtoE1/.ruleQtoE2/.ER->53.14*10^-6//Simplify}]
  ,{p,1,3}]
,{i,1,Nops}]
*)

(* Exit with success exit code *)
Exit[0];
