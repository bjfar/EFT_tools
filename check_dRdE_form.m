#!/usr/local/bin/MathKernel -script

(* Compute recoil spectra for a specified isotope of Xenon *)

(*Args:*)
(*1 - Path to "v6/dmformfactor-V6.m"*)
(*2 - Xenon isotope (e.g. "131")*)
(*3 - *)


args = $CommandLine[[4;;]]
Print["Running with command line arguments: "<>ToString[args]]
path = args[[1]]
Print["Setting working directory to: "<>path]
isotope = ToExpression[args[[2]]]
SetDirectory[path]
<<"v6/dmformfactor-V6.m";

(* Model Setup *)
SetJChi[1/2](* WIMP spin *)
SetMChi[MWIMP GeV] (* WIMP mass in GeV *)
densmatrix="default";  (* nuclear density matrix (using default density matrices) *)
bFM="default"; (* oscillation parameter (default: approximate formulae used) *)
IsotopeA=isotope;
Print["Setting A="<>ToString[IsotopeA]]
SetIsotope[54,IsotopeA,bFM,densmatrix] (*[Z,A,bFM,filename]*)
ZeroCoeffs[];(*Reset all coefficients*)
SetCoeffsNonrel[1,c1p,"p"]
SetCoeffsNonrel[1,c1n,"n"]
SetCoeffsNonrel[2,c2p,"p"]
SetCoeffsNonrel[2,c2n,"n"]
SetCoeffsNonrel[3,c3p,"p"]
SetCoeffsNonrel[3,c3n,"n"]
SetCoeffsNonrel[4,c4p,"p"]
SetCoeffsNonrel[4,c4n,"n"]
SetCoeffsNonrel[5,c5p,"p"]
SetCoeffsNonrel[5,c5n,"n"]
SetCoeffsNonrel[6,c6p,"p"]
SetCoeffsNonrel[6,c6n,"n"]
SetCoeffsNonrel[7,c7p,"p"]
SetCoeffsNonrel[7,c7n,"n"]
SetCoeffsNonrel[8,c8p,"p"]
SetCoeffsNonrel[8,c8n,"n"]
SetCoeffsNonrel[9,c9p,"p"]
SetCoeffsNonrel[9,c9n,"n"]
SetCoeffsNonrel[10,c10p,"p"]
SetCoeffsNonrel[10,c10n,"n"]
SetCoeffsNonrel[11,c11p,"p"]
SetCoeffsNonrel[11,c11n,"n"]
SetCoeffsNonrel[12,c12p,"p"]
SetCoeffsNonrel[12,c12n,"n"]
SetCoeffsNonrel[13,c13p,"p"]
SetCoeffsNonrel[13,c13n,"n"]
SetCoeffsNonrel[14,c14p,"p"]
SetCoeffsNonrel[14,c14n,"n"]
SetCoeffsNonrel[15,c15p,"p"]
SetCoeffsNonrel[15,c15n,"n"]
(* Set non-zero NR EFT coefficient values *)


(* Parameters *)
mNucleon=0.938 GeV;
NT=1/(IsotopeA* mNucleon);(*Number of target nuclei per detector mass*)
Centimeter=(10^13 Femtometer);
rhoDM=0.3GeV/Centimeter^3;(*Local dark matter density*)
ve=232 KilometerPerSecond;(*Earth's velocity in galactic rest frame*)
v0=220 KilometerPerSecond;(*Mean WIMP speed in galactic rest frame*)
vesc=550 KilometerPerSecond;
SetHalo["MBcutoff"];

(* Event rate in symbolic form *)
dRdE=EventRate[NT,rhoDM,qGeV,ve,v0,vesc]; (* as function of q, not ER *)


(* Produce recoil spectra for several WIMP masses *)
m\[Chi]={5,10,50,100,500};(*GeV*)
mT = IsotopeA*mNucleon;
\[Mu]T=Table[(m\[Chi][[i]] GeV mT)/(m\[Chi][[i]] GeV+mT),{i,1,Length[m\[Chi]]}];
ERmax =Table[ 2(((\[Mu]T[[i]]^2) (v^2) )/mT)/GeV/.v->ve,{i,1,Length[m\[Chi]]}];(*just using earth speed to get rough number, and get rid of GeV units for plotting*)
ERmaxesc = Table[2(((\[Mu]T[[i]]^2) (v^2) )/mT)/GeV/.v->vesc,{i,1,Length[m\[Chi]]}];
(*generalfunc= Table[2500KilogramDay*fEdRdE[ER,m\[Chi]\[LeftDoubleBracket]i\[RightDoubleBracket],mT/GeV],{i,1,Length[m\[Chi]]}];*)


(* General function to compute recoil spectra (couplings not yet replaced with values) *)
generalfunc= Table[7800 KilogramDay*dRdE GeV/.qGeV->Sqrt[2mTv ER]/.MWIMP->m\[Chi][[i]]/.mTv->mT/GeV,{i,1,Length[m\[Chi]]}];

curvenames=Table["m\[Chi]="<>ToString[m\[Chi][[i]]],{i,1,Length[m\[Chi]]}]

(* Replacement rules for setting the couplings one by one *)
coefflist={c1p,c1n,c2p,c2n,c3p,c3n,c4p,c4n,c5p,c5n,c6p,c6n,c7p,c7n,c8p,c8n,c9p,c9n,c10p,c10n,c11p,c11n,c12p,c12n,c13p,c13n,c14p,c14n,c15p,c15n};
resttozero={c1p->0,c1n->0,c2p->0,c2n->0,c3p->0,c3n->0,c4p->0,c4n->0,c5p->0,c5n->0,c6p->0,c6n->0,c7p->0,c7n->0,c8p->0,c8n->0,c9p->0,c9n->0,c10p->0,c10n->0,c11p->0,c11n->0,c12p->0,c12n->0,c13p->0,c13n->0,c14p->0,c14n->0,c15p->0,c15n->0};
Nops=15

(* Full list of all the rules we want *)
rules=Table[{{coefflist[[2*i-1]]->1},{coefflist[[2*i]]->1},{coefflist[[2*i-1]]->1,coefflist[[2*i]]->1}},{i,1,Nops}]
Monitor[funcsinter=Table[Table[Table[generalfunc[[mWIMP]]/.rules[[i]][[p]]/.resttozero,{mWIMP,1,Length[m\[Chi]]}],{p,1,3}],{i,1,Nops}],{i,p,mWIMP}];

(* plot titles *)
titles=Table[{ToString[coefflist[[2*i-1]]],ToString[coefflist[[2*i]]],ToString[coefflist[[2*i-1]]]<>"="<>ToString[coefflist[[2*i]]]},{i,1,Nops}]

(* Apply all the rules so that only the recoil energies are left unspecified*)
Monitor[funcs=Table[Table[Table[With[{i=i,p=p,mWIMP=mWIMP},(funcsinter[[i]][[p]][[mWIMP]]/.ER->#)&],{mWIMP,1,Length[m\[Chi]]}],{p,1,3}],{i,1,Nops}],{i,p,mWIMP}];

(* Specifiy recoil energies at which to evaluate the function *)
Earr=Table[N[ER*10^-6],{ER,0 ,1000,1}];

(* Apply the analytic expressions to the chosen recoil energies and obtain final numerical tables *)
Monitor[plotdata=Table[Table[Transpose[Join[{Earr*10^6},Table[10^-6*funcs[[i]][[p]][[m]]/@Earr,{m,1,Length[m\[Chi]]}]]],{p,1,3}],{i,1,15}],{i,p,m}];

(* Save data tables to file *)
Do[Do[Export["/home/farmer/mathematica/DMFormFactor_13086288/EFTcoeffplotdata/Xe"<>ToString[IsotopeA]<>"/Xe"<>ToString[IsotopeA]<>"_"<>titles[[i]][[p]]<>".dat",plotdata[[i]][[p]],"Table"],{i,1,Nops}],{p,1,3}];


(* Report on which operators vanish for this isotope *)
Do[Do[Print[{"Vanishes? ",0==dRdE GeV/.qGeV->Sqrt[2mTv ER]/.MWIMP->500/.mTv->mT/GeV/.rules[[i]][[p]]/.resttozero/.ER->53.14*10^-6,rules[[i]][[p]]}],{p,1,3}],{i,1,Nops}]



