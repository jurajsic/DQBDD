#!/bin/sh

benchmarkFiles="benchmarks/dqbf18/*"

timeLength="10m"

for pathToFile in $benchmarkFiles; do
    fileName=$(basename $pathToFile)
    echo $fileName
    lineFormat=$fileName",%e,%x"
    time -q -o "results/simplesolver.csv" -a --format=$lineFormat timeout $timeLength solvers/solver 0 $pathToFile
    time -q -o "results/treesolver.csv" -a --format=$lineFormat timeout $timeLength solvers/solver 1 $pathToFile
    time -q -o "results/hqsold.csv" -a --format=$lineFormat timeout $timeLength solvers/hqs $pathToFile
    time -q -o "results/hqslq.csv" -a --format=$lineFormat timeout $timeLength solvers/HQSnp $pathToFile
done