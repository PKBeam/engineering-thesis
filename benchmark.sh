#!/bin/bash
echo N,S,T1,T2a,T2b,T3,T4,E1,E2a,E2b,E3,E4 > benchmark.csv
for s in 1 0.5 0.1 0.05 0.01 0.005
do
    for N in 150 200
    #for N in {1..100}
    do
	echo benchmarking N = $N, step size $s...
        ./Thesis\ B.wls $N $s >> benchmark.csv
    done
done
