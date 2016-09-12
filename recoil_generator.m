#!/usr/local/bin/MathKernel -script

(* Compute recoil spectra for a specified isotope of Xenon *)

(*Args:*)
(*1 - Absolute path to "v6/dmformfactor-V6.m"*)
(*2 - Relative output path*)
(*3 - Xenon isotope (e.g. "131")*)
(*4 - flag to use relativistic operator coefficients*)
(*5 - flag to switch form factor treatment to traditional Helm form factors rather than derived from nuclear density matrices *) 

(* Make output more readable in the terminal *)
(* SetOptions["stdout", FormatType->InputForm] *)

args = Join[$CommandLine[[4;;]],Table[{},{i,1,10}]];
(* Print["Running with command line arguments: "<>ToString[args]]; *)
outpath = Directory[]<>"/"<>args[[2]];
Print["Setting output directory to: "<>outpath];
If[!DirectoryQ[outpath],
  Print["\nOutput directory "<>outpath<>" does not exist! Please create it and try again.\n"];
  Exit[1];
];

path = args[[1]];
Print["Looking for DMFormFactor in directory: "<>path]
If[!DirectoryQ[path],
  Print["\nSupplied DMFormFactor directory "<>path<>" does not exist! Please check it for typos and try again\n"];
  Exit[1];
];
SetDirectory[path];
(* SetDirectory["/home/farmer/mathematica/DMFormFactor_13086288"] (* for testing in notebook *)*)
Check[
  <<"v6/dmformfactor-V6.m";
 ,
  Print["...Error loading DMFormFactor package. Please check the path (and perhaps version; this script requires v6) and try again"];
  Exit[1];
]
Print["...Found."];
isotope = ToExpression[args[[3]]];
Otype = args[[4]]
userel = SameQ[Otype,"R"]; (* Switch to using coefficients of relativistic operators (spin 1/2 WIMP only) *)
checkrel = SameQ[args[[5]],"check"]; (* Compute coefficients of relativistic operators using their non-relativistic reductions (to check that I understand the mapping correctly *)
helm = SameQ[args[[6]],"UseHelm"]; (* Use "traditional" Helm form factors rather than ones derived from nuclear density matrices *)


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
m\[Chi]={5,10,50}; (* GeV *)
(*m\[Chi]={3,4,5,6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,30,35,40,45,50,60,70,80,100,200,300,500,700,800,1000,2000}; (* GeV *)*)
(*m\[Chi]={3,5,6,7,8,9,10,15,20,30,50,100,300,500,700,800,1000,2000,3000,5000}; (* GeV *) *)
(*m\[Chi]={5,50,500,5000}; (* GeV *) *)
spin=1/2 (*leave hardcoded for now*)

(* Exposure to use for rate calculations. Can be easily scaled to something different in later analysis steps. *) 
benchmarkExposure = 7800 KilogramDay;

(* Turn on handler to abort evaluation if anything weird happens *)
messageHandler = If[Last[#], Print["An error has occurred in 'recoil_generator.m'; aborting evaluation..."] Exit[1];]&
Internal`AddHandler["Message", messageHandler]

Print["**************************************"    ]
Print[" Generating recoil spectra for Xe"<>ToString[isotope]<>"..."]
Print[" ...for spin "<>ToString[spin,InputForm]<>" WIMP"    ]
Print[" ...for masses (GeV) "<>ToString[m\[Chi]]  ]
Print[" ...for operator type "<>Otype             ]
If[helm,
Print[" ...with traditional Helm form factors"    ]]
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
     SetCoeffsRel[1,c1p,"p"];
     SetCoeffsRel[1,c1n,"n"];
     SetCoeffsRel[2,c2p,"p"];
     SetCoeffsRel[2,c2n,"n"];
     SetCoeffsRel[3,c3p,"p"];
     SetCoeffsRel[3,c3n,"n"];
     SetCoeffsRel[4,c4p,"p"];
     SetCoeffsRel[4,c4n,"n"];
     SetCoeffsRel[5,c5p,"p"];
     SetCoeffsRel[5,c5n,"n"];
     SetCoeffsRel[6,c6p,"p"];
     SetCoeffsRel[6,c6n,"n"];
     SetCoeffsRel[7,c7p,"p"];
     SetCoeffsRel[7,c7n,"n"];
     SetCoeffsRel[8,c8p,"p"];
     SetCoeffsRel[8,c8n,"n"];
     SetCoeffsRel[9,c9p,"p"];
     SetCoeffsRel[9,c9n,"n"];
     SetCoeffsRel[10,c10p,"p"];
     SetCoeffsRel[10,c10n,"n"];
     SetCoeffsRel[11,c11p,"p"];
     SetCoeffsRel[11,c11n,"n"];
     SetCoeffsRel[12,c12p,"p"];
     SetCoeffsRel[12,c12n,"n"];
     SetCoeffsRel[13,c13p,"p"];
     SetCoeffsRel[13,c13n,"n"];
     SetCoeffsRel[14,c14p,"p"];
     SetCoeffsRel[14,c14n,"n"];
     SetCoeffsRel[15,c15p,"p"];
     SetCoeffsRel[15,c15n,"n"];
     SetCoeffsRel[16,c16p,"p"];
     SetCoeffsRel[16,c16n,"n"];
     SetCoeffsRel[17,c17p,"p"];
     SetCoeffsRel[17,c17n,"n"];
     SetCoeffsRel[18,c18p,"p"];
     SetCoeffsRel[18,c18n,"n"];
     SetCoeffsRel[19,c19p,"p"];
     SetCoeffsRel[19,c19n,"n"];
     SetCoeffsRel[20,c20p,"p"];
     SetCoeffsRel[20,c20n,"n"];
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
     SetCoeffsNonrel[1,c1p,"p"];
     SetCoeffsNonrel[1,c1n,"n"];
     SetCoeffsNonrel[2,c2p,"p"];
     SetCoeffsNonrel[2,c2n,"n"];
     SetCoeffsNonrel[3,c3p,"p"];
     SetCoeffsNonrel[3,c3n,"n"];
     SetCoeffsNonrel[4,c4p,"p"];
     SetCoeffsNonrel[4,c4n,"n"];
     SetCoeffsNonrel[5,c5p,"p"];
     SetCoeffsNonrel[5,c5n,"n"];
     SetCoeffsNonrel[6,c6p,"p"];
     SetCoeffsNonrel[6,c6n,"n"];
     SetCoeffsNonrel[7,c7p,"p"];
     SetCoeffsNonrel[7,c7n,"n"];
     SetCoeffsNonrel[8,c8p,"p"];
     SetCoeffsNonrel[8,c8n,"n"];
     SetCoeffsNonrel[9,c9p,"p"];
     SetCoeffsNonrel[9,c9n,"n"];
     SetCoeffsNonrel[10,c10p,"p"];
     SetCoeffsNonrel[10,c10n,"n"];
     SetCoeffsNonrel[11,c11p,"p"];
     SetCoeffsNonrel[11,c11n,"n"];
     SetCoeffsNonrel[12,c12p,"p"];
     SetCoeffsNonrel[12,c12n,"n"];
     SetCoeffsNonrel[13,c13p,"p"];
     SetCoeffsNonrel[13,c13n,"n"];
     SetCoeffsNonrel[14,c14p,"p"];
     SetCoeffsNonrel[14,c14n,"n"];
     SetCoeffsNonrel[15,c15p,"p"];
     SetCoeffsNonrel[15,c15n,"n"];
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
Print["Computing symbolic event rate..."];
DoSilent[
  dRdE=EventRate[NT,rhoDM,qGeV,ve,v0,vesc]; (* as function of q, not ER *)
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
generalfunc= benchmarkExposure*dRdE GeV /. ruleQtoE1 /. ruleQtoE2;

curvenames=Table["m\[Chi]="<>ToString[m\[Chi][[i]]],{i,1,Length[m\[Chi]]}];

(* Replacement rules for setting the couplings one by one *)
coefflist={c1p,c1n,c2p,c2n,c3p,c3n,c4p,c4n,c5p,c5n,c6p,c6n,c7p,c7n,c8p,c8n,c9p,c9n,c10p,c10n,c11p,c11n,c12p,c12n,c13p,c13n,c14p,c14n,c15p,c15n,c16p,c16n,c17p,c17n,c18p,c18n,c19p,c19n,c20p,c20n};

coefflistnoiso={c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c14,c15,c16,c17,c18,c19,c20};

resttozero={c1p->0,c1n->0,c2p->0,c2n->0,c3p->0,c3n->0,c4p->0,c4n->0,c5p->0,c5n->0,c6p->0,c6n->0,c7p->0,c7n->0,c8p->0,c8n->0,c9p->0,c9n->0,c10p->0,c10n->0,c11p->0,c11n->0,c12p->0,c12n->0,c13p->0,c13n->0,c14p->0,c14n->0,c15p->0,c15n->0,c16p->0,c16n->0,c17p->0,c17n->0,c18p->0,c18n->0,c19p->0,c19n->0,c20p->0,c20n->0,
c1->0,c2->0,c3->0,c4->0,c5->0,c6->0,c7->0,c8->0,c9->0,c10->0,c11->0,c12->0,c13->0,c14->0,c15->0,c16->0,c17->0,c18->0,c19->0,c20->0};

(* Replacement rules for recreating relativistic operators by setting their non-relativistic reductions correctly (see table 1 of documentation-standalone.pdf) *)
(* Result should be equivalent to setting the various d's to 1 *)

(* Note: "The coefficients d_j are dimensionless, by inserting appropriate powers of the user-defined scale m_M (which usually would be known from the context of the theory), set by the user function SetMM. This scale is set by default to be m M = m V ≡ 246.2 GeV" *)
(* So following the above, I will use the default mM=246.2 for now *)
mM = 246.2 GeV;
mX = MWIMP GeV;
mN = mNucleon;
reductionrules=
{
  (*1*) {c1->1},
  (*2*) {c10->1},
  (*3*) {c11->-mN/mX},
  (*4*) {c6-> -mN/mX},
  (*5*) {c1->1},
  (*6*) {c1-> qGeV^2/(2*mN*mM), 
         c3->-2*mN/mM, 
         c4-> 2*mN^2/(mM*mX) * qGeV^2/mN^2,  
         c6->-2*mN^2/(mM*mX)},
  (*7*) {c7->-2,
         c9->2*mN/mX},
  (*8*) {c10->2*mN/mM},
  (*9*) {c1->-qGeV^2/(2*mX*mM),
         c5->2*mN/mM,
         c4->-2*mN/mM * qGeV^2/mN^2,
         c6-> 2*mN/mM},
  (*10*) {c4->4*qGeV^2/mM^2,
          c6->-mN^2/mM^2},
  (*11*) {c9->4*mN/mM},
  (*12*) {c10->-mN/mX * qGeV^2/mM^2,
          c12->-4*qGeV^2/mM^2,
          c15->-4*mN^2/mM^2},
  (*13*) {c8->2,
          c9->2},
  (*14*) {c9->-4*mN/mM},
  (*15*) {c4->-4},
  (*16*) {c13->4*mN/mM},
  (*17*) {c11->2*mN/mM},
  (*18*) {c11->qGeV^2/mM^2,
          c15->4*mN^2/mM^2},
  (*19*) {c14->-4*mN/mM},
  (*20*) {c6->4*mN^2/mM^2}
};


(* Full list of all the rules we want *)
If[userel && checkrel,
  Print["Applying coefficient replacement rules (using manual non-relativistic reductions)..."];
  ponly=Table[coefflist[[2*i-1]]->coefflistnoiso[[i]],{i,1,15}];
  nonly=Table[coefflist[[2*i]]->coefflistnoiso[[i]],{i,1,15}];
  nandp=Join[ponly,nonly];
  (* rules are: first, remove p and/or n from names of non-zero coefficients (rest go to zero). Second, apply one of the non-relativistic reduction rules *)
  rules=Table[ { Join[ponly,reductionrules[[i]]],
                 Join[nonly,reductionrules[[i]]],
                 Join[nandp,reductionrules[[i]]] } ,{i,1,Nops}];
  , (* else do normal thing *)
  Print["Applying coefficient replacement rules..."];
  rules=Table[{{coefflist[[2*i-1]]->1},{coefflist[[2*i]]->1},{coefflist[[2*i-1]]->1,coefflist[[2*i]]->1}},{i,1,Nops}];
]

funcsinter=Table[Table[ 
  WriteString["stdout", "  Evaluating WIMP mass= "<>ToString[m\[Chi][[m]]]<>", Operator="<>ToString[i]<>"          ",  "\r"];
  Table[
    temp=(generalfunc//.rules[[i]][[p]]/.resttozero/.MWIMP->m\[Chi][[m]]);
    (* Simplify the result if all the GeV's haven't already cancelled out *)
    If[ FreeQ[temp, GeV],
      temp,
      Simplify[temp]
    ]
  ,{p,1,3}]
,{m,1,Length[m\[Chi]]}],{i,1,Nops}];
Print["  Complete!                                "];

(* plot titles *)
titles=Table[{ToString[coefflist[[2*i-1]]],ToString[coefflist[[2*i]]],ToString[coefflist[[2*i-1]]]<>"="<>ToString[coefflist[[2*i]]]},{i,1,Nops}];

(* Convert expressions to functions so that recoil energies can be supplied as arguments*)
Print["Converting expressions to functions..."]
funcs=Table[Table[
  WriteString["stdout", "  Evaluating WIMP mass= "<>ToString[m\[Chi][[m]]]<>", Operator="<>ToString[i]<>"          ", "\r"];
  Table[
    With[{i=i,p=p,m=m},(funcsinter[[i]][[m]][[p]]/.ER->#)&]
  ,{p,1,3}]
,{m,1,Length[m\[Chi]]}],{i,1,Nops}];
Print["  Complete!                                "];

(* Specifiy recoil energies at which to evaluate the function *)
Earr=Table[N[ER*10^-6],{ER,1,1000,1}];

(* Apply the analytic expressions to the chosen recoil energies and obtain final numerical tables *)
Print["Computing numerical results at "<>ToString[Length[Earr]]<>" recoil energies between "<>ToString[FortranForm[First[Earr]]]<>" and "<>ToString[FortranForm[Last[Earr]]]<>" GeV"];
plotdata=
 Table[
   WriteString["stdout", "  Evaluating operator="<>ToString[i]<>"          ", "\r"];
   Table[
      Transpose[Join[{Earr*10^6},
                  Table[
                        10^-6*funcs[[i]][[m]][[p]]/@Earr
                       ,{m,1,Length[m\[Chi]]}
                       ]
                    ]
               ]
   ,{p,1,3}]
 ,{i,1,Nops}];
Print["  Complete!                                "];

(* Save data tables to file *)
extendedoutpath=outpath<>"/Xe"<>ToString[IsotopeA]
If[!DirectoryQ[extendedoutpath], CreateDirectory[extendedoutpath]]
Do[
  Do[
    filename="Xe"<>ToString[IsotopeA]<>"_"<>Tag<>"_"<>titles[[i]][[p]]<>".dat";
    Export[extendedoutpath<>"/"<>filename, plotdata[[i]][[p]], "Table"]
  ,{i,1,Nops}]
,{p,1,3}];

(* Report on which operators vanish for this isotope *)
Do[
  Do[
    Print[{"Vanishes? ",0==dRdE GeV//.rules[[i]][[p]]/.resttozero/.MWIMP->500/.ruleQtoE1/.ruleQtoE2/.ER->53.14*10^-6//Simplify}]
  ,{p,1,3}]
,{i,1,Nops}]

(* Exit with success exit code *)
Exit[0];
