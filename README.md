## Applications of Local Signal Approximations  to Differential Equations

The chromatic derivative is a linear differential operator originating from the field of signal analysis. The aim of this thesis is to investigate applications of the chromatic derivative to the numerical solution of differential equations. 

We propose methods to solve linear ordinary differential equations using chromatic expansions, which are analogues of Taylor expansions for the chromatic derivative. This has many potential applications, such as in the implementation of differential equation solvers in mathematical software.

Finally, we implement these methods in software using the Wolfram language and evaluate them against known numerical differential equation solvers, focusing on execution time and accuracy. 

### Usage

#### Interactive
Open the Mathematica notebook `Thesis B.nb` and run it.  
The differential equation problem can be modified inside the notebook.

#### Manual
Edit `Thesis B.wls` (resp. `Thesis B control.wls`) to set up the differential equation problem.  
Run `benchmark.sh` (resp. `benchmark control.sh`). Benchmark output will be in `benchmark.csv`.
