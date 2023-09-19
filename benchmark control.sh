#!/bin/bash
echo A,T,E > benchmark.csv
for N in {1..20}
do
    ./Thesis\ B\ Control.wls $N >> benchmark.csv
done
