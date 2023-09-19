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
Off[NIntegrate::ncvb]
Off[NIntegrate::slwcon]
Off[DecimalForm::reqsigz]
Off[N::meprec]
Off[General::munfl]

(* ::Package:: *)

(* ::Input::Initialization:: *)
(* Thesis Project: 
Approximating solutions to differential equations using chromatic expansions.

By Vincent Liu (z5258281).
Supervisor: Aleksandar Ignjatovi\[CAcute]. 

In this notebook, we will solve a differential equation by considering the chromatic expansion of the solution. This method will be analysed,    benchmarked and compared to other methods.
*)

(* ::Input::Initialization:: *)
(* Parameters for the chromatic approximation process.

Nd - the precision of the expansion (i.e. the series truncates after Nd terms).
pt - the set of points over which the solution is approximated on.
K - the number of points in the set pt. 
*)
Nd = ToExpression[$ScriptCommandLine[[2]]];
step = ToExpression[$ScriptCommandLine[[3]]];
pt = Table[i, {i, 0, 20, step}];
K = Length[pt];
(*Print["N = ", Nd, ", pt: (", pt[[1]], ", ", pt[[K]], "), step = ", (pt[[K]]-pt[[1]])/(K-1), ", K = ", K]*)

(* ::Input::Initialization:: *)
(* Recurrence co-efficients for Legendre polynomials. *)
r[-1]:=0; r[-2]:=0; r[n_]=( \[Pi](n+1))/Sqrt[4(n+1)^2-1];
(* Basis functions for the chromatic approximations. *)
J[n_][t_] = Re[Sqrt[2n+1]SphericalBesselJ[n, \[Pi] t]];

(* The m-th derivative of J[n,t] for t = 0 to p. 

These are calculated using the recurrence coefficients defined above.
If we try to differentiate J with Mathematica we encounter numerically unstable matrices later. 
*)
DJ[0][n_,t_]:=J[n][t];
DJ[m_][n_,t_]:= r[n-1] DJ[m-1][n-1, t]- r[n] DJ[m-1][n+1, t];


(* ::Input::Initialization:: *)
(* The differential equation we'll try to solve is
     L[x, t] = f[t]. 

The parameter p represents the order of the differential equation.
*)
p = 3;
\[Omega] = 2;
f[t_]=t^2;
(*L[x_, t_] = x''[t]-t x[t];*)
(*L[x_, t_] = t^2 x''[t] + t x'[t] + (t^2 - 1) x[t];*)
L[x_, t_] =t^3x'''[t]-t^2x''[t]-t x'[t]+x[t];
(*L[x_, t_] =x''[t]+\[Omega]^2 x[t];*)

(* The differential operator L evaluated on J[n,t].

We need to specify this separately because we need Mathematica to differentiate J using our function DJ defined above. 
*)
(*LJ[n_][t_] = DJ[2][n, t] - t J[n][t];*)
(*LJ[n_][t_] = t^2 DJ[2][n, t] + t DJ[1][n, t] + (t^2 - 1) J[n][t];*)
LJ[n_][t_] = t^3DJ[3][n,t]-t^2DJ[2][n,t]-t DJ[1][n,t]+J[n][t];
(*LJ[n_][t_] = DJ[2][n,t]+\[Omega]^2 J[n][t];*)

(* Initial conditions for the differential equation.
They are provided in the form
	(x^(m))[T[[m]]] = c[[m]]. 
*)
T= Table[i, {i, 0, p-1}];
c =Table[i, {i, 0, p-1}];
T[[1]]=0;
c[[1]]=0;
T[[2]]=1;
c[[2]]=1;
T[[3]]=1;
c[[3]]=1;
Dx[m_][t_]=D[x[t], {t,m}];
ICs = Table[Dx[i-1][T[[i]]]== c[[i]], {i, 1, p}];


(* ::Input::Initialization:: *)
(* Set up matrices for linear algebra solutions. 
We want to find the closest solution x to Ax = b. 
The vector x represents the approximated coefficients of the solution's chromatic expansion. 
*)
Tms = RepeatedTiming[
A = ConstantArray[0, {p+K, Nd+1}];
Quiet[Do[A[[i, j]]=DJ[i-K-1][j-1,T[[i-K]]],{i, K+1, K+p}, {j,1,Nd+1}]];
Do[A[[i,j]]=LJ[j-1][pt[[i]]], {i, 1, K}, {j,1,Nd+1}];

b = ConstantArray[0, {p+K}];
Do[b[[i]]=c[[i-K]], {i, K+1, K+p}];
Do[b[[i]]=f[pt[[i]]], {i, 1, K}];

nA= N[A];nb = N[b];
][[1]];


(* ::Input::Initialization:: *)
(* 0. Find the exact solution *)
S = DSolve[Join[{L[x, t] == f[t]}, ICs],  x, {t, pt[[1]], pt[[K]]}];


(* ::Input::Initialization:: *)
(* 1. Find an approximate solution with NDSolve *)
TS1=RepeatedTiming[
S1 = NDSolve[Join[{L[x, t] == f[t]}, ICs],  x, {t, pt[[1]], pt[[K]]}];
][[1]];


(* ::Input::Initialization:: *)
(* 2a. Manually solve with Least Squares *)
TS2a=RepeatedTiming[
nAT = Transpose[nA];
nATA = nAT . nA;
nATb = nAT . nb;
R = 10^-5 IdentityMatrix[Length[nATA]];
X2a = Inverse[nATA+R] . nATb;
S2a[t_] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\(X2a[\([n + 1]\)] \(J[n]\)[t]\)\);
][[1]];


(* ::Input::Initialization:: *)
(* 2b. Manually solve with Least Norm *)
TS2b=RepeatedTiming[
nAT = Transpose[nA];
nAAT = nA . nAT;
R = 10^-5 IdentityMatrix[Length[nAAT]];
X2b = nAT . Inverse[nAAT+R] . nb;
S2b[t_] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\(X2b[\([n + 1]\)] \(J[n]\)[t]\)\);
][[1]];


(* ::Input::Initialization:: *)
(* 3. Solve with Mathematica's LeastSquares *)
TS3=RepeatedTiming[
X3=LeastSquares[nA,nb];
S3[t_] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\(X3[\([n + 1]\)] \(J[n]\)[t]\)\);
][[1]];


(* ::Input::Initialization:: *)
(* 4. Solve with NMinimize *)
TS4=RepeatedTiming[
XMs=Table[XM[n],  {n, 0, Nd}];
Quiet[M4=NMinimize[ \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 0\), \(p - 1\)]
\*SuperscriptBox[\((
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\((XM[n] \(DJ[i]\)[n, T[\([\)\(i + 1\)\(]\)]])\) - c[\([\)\(i + 1\)\(]\)])\), \(2\)]\)+\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(K\)]
\*SuperscriptBox[\((
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\((XM[n] \(LJ[n]\)[pt[\([\)\(i\)\(]\)]])\) - f[pt[\([\)\(i\)\(]\)]])\), \(2\)]\)+10^-4 \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\((
\*SuperscriptBox[\(XM[n]\), \(2\)])\)\),XMs];
];
X4=Re[XMs/.Last[M4]];
S4[t_] = \!\(
\*UnderoverscriptBox[\(\[Sum]\), \(n = 0\), \(Nd\)]\(X4[\([\)\(n + 1\)\(]\)] \(J[n]\)[t]\)\);
][[1]];


(* ::Input::Initialization:: *)
(* Plot solutions *)
(*Plot[Evaluate[x[t]/. S],{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[Evaluate[x[t]/. S1],{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[S2[t],{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[S3[t],{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[S4[t],{t,pt[[1]],pt[[K]]},PlotRange->All]*)


(* ::Input::Initialization:: *)
(* Plot the error *)
E1=Evaluate[x[t]/. S] - Evaluate[x[t]/. S1];
E2a=Evaluate[x[t]/. S] - S2a[t];
E2b=Evaluate[x[t]/. S] - S2b[t];
E3=Evaluate[x[t]/. S] - S3[t];
E4=Evaluate[x[t]/. S] - S4[t];

(*Plot[E1,{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[E2,{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[E3,{t,pt[[1]],pt[[K]]},PlotRange->All]
Plot[E4,{t,pt[[1]],pt[[K]]},PlotRange->All]*)


(* ::Input::Initialization:: *)
(* Evaluate the error (as L2 norm) *)
ES=Sqrt[NIntegrate[Evaluate[x[t]/. S], {t,pt[[1]],pt[[K]]}]][[1]];
ES1=Sqrt[NIntegrate[E1^2, {t,pt[[1]],pt[[K]]}]][[1]];
ES2a=Sqrt[NIntegrate[E2a^2, {t,pt[[1]],pt[[K]]}]][[1]];
ES2b=Sqrt[NIntegrate[E2b^2, {t,pt[[1]],pt[[K]]}]][[1]];
ES3=Sqrt[NIntegrate[E3^2, {t,pt[[1]],pt[[K]]}]][[1]];
ES4=Sqrt[NIntegrate[E4^2, {t,pt[[1]],pt[[K]]}]][[1]];


(* ::Input::Initialization:: *)
(* Print out the evaluation results. 
Note: the chromatic solutions need to add the time for matrix setup.
*)
Print[Nd,",",step,",",1000*TS1,",",1000*(Tms+TS2a),",",1000*(Tms+TS2b),",",1000*(Tms+TS3),",",1000*TS4,",",DecimalForm[ES1],",",DecimalForm[ES2a],",",DecimalForm[ES2b],",",DecimalForm[ES3],",",DecimalForm[ES4]]



