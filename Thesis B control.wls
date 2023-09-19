#!/Applications/Mathematica.app/Contents/MacOS/wolframscript
SetOptions[$Output, FormatType -> OutputForm]
Off[Power::infy]
Off[Infinity::indet]
Off[General::stop]
Off[NDSolve::nlnum]
Off[NDSolve::ndnum]
Off[NDSolve::berr]
Off[NDSolve::ndsz]
Off[Divide::infy]
Off[Inverse::luc]
Off[NIntegrate::inumr]
Off[NIntegrate::slwcon]
Off[NIntegrate::ncvb]

A = ToExpression[$ScriptCommandLine[[2]]];
pt = Table[i, {i, -10, 10, 0.5} ];
K = Length[pt];

(* ::Input::Initialization:: *)
(* The differential equation we'll try to solve is
     L[x, t] = f[t]. 

The parameter p represents the order of the differential equation.
*)
p = 2;
\[Omega] = 2;
f[t_]=0;
L[x_, t_] = x''[t]-t x[t];
(*L[x_, t_] = t^2 x''[t] + t x'[t] + (t^2 - 1) x[t];*)
(*L[x_, t_] =t^3x'''[t]-t^2x''[t]-t x'[t]+x[t];*)
(*L[x_, t_] =x''[t]+\[Omega]^2 x[t];*)

(* The differential operator L evaluated on J[n,t].

We need to specify this separately because we need Mathematica to differentiate J using our function DJ defined above. 
*)
LJ[n_][t_] = DJ[2][n, t] - t J[n][t];
(*LJ[n_][t_] = t^2 DJ[2][n, t] + t DJ[1][n, t] + (t^2 - 1) J[n][t];*)
(*LJ[n_][t_] = t^3DJ[3][n,t]-t^2DJ[2][n,t]-t DJ[1][n,t]+J[n][t];*)
(*LJ[n_][t_] = DJ[2][n,t]+\[Omega]^2 J[n][t];*)

(* Initial conditions for the differential equation.
They are provided in the form
	(x^(m))[T[[m]]] = c[[m]]. 
*)
T= Table[i, {i, 0, p-1}];
c =Table[i, {i, 0, p-1}];
T[[1]]=0;
c[[1]]=0;
T[[2]]=0;
c[[2]]=1;
(*T[[3]]=1;
c[[3]]=1;*)
Dx[m_][t_]=D[x[t], {t,m}];
ICs = Table[Dx[i-1][T[[i]]]== c[[i]], {i, 1, p}];

(* ::Input::Initialization:: *)
(* 0. Find the exact solution *)
S = DSolve[Join[{L[x, t] == f[t]}, ICs],  x, {t, pt[[1]], pt[[K]]}];

(* ::Input::Initialization:: *)
(* 1. Find an approximate solution with NDSolve *)
TS1=RepeatedTiming[
S1 = NDSolve[Join[{L[x, t] == f[t]}, ICs],  x, {t, pt[[1]], pt[[K]]}, AccuracyGoal -> A];
][[1]];

E1=Evaluate[x[t]/. S] - Evaluate[x[t]/. S1];
ES1=Sqrt[NIntegrate[E1^2, {t,pt[[1]],pt[[K]]}]][[1]];

(* ::Input::Initialization:: *)
(* Print out the evaluation results. 
Note: the chromatic solutions need to add the time for matrix setup.
*)
Print[A,",",1000*TS1,",",DecimalForm[ES1]]



