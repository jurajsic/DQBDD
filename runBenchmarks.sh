#!/bin/sh

benchmarkFiles="benchmarks/dqbf18/*"

timeLength="10m"

mkdir results
mkdir results/simpleOutput
mkdir results/treeOutput
mkdir results/hqsoldOutput
mkdir results/hqslqOutput

for pathToFile in $benchmarkFiles; do
    fileName=$(basename $pathToFile)
    echo "Running solvers for $fileName"
    lineFormat=$fileName",%e,%x"
    echo  "Simple BDD solver"
    time -q -o "results/simplesolver.csv" -a --format=$lineFormat timeout $timeLength solvers/solver 0 $pathToFile > "results/simpleOutput/$fileName"
    echo "Tree BDD solver"
    time -q -o "results/treesolver.csv" -a --format=$lineFormat timeout $timeLength solvers/solver 1 $pathToFile > "results/treeOutput/$fileName"
    echo "Old HQS with older HQSpre"
    time -q -o "results/hqsold.csv" -a --format=$lineFormat timeout $timeLength solvers/hqs $pathToFile > "results/hqsoldOutput/$fileName"
    echo "HQS with quantifier localisation and newer HQSpre"
    time -q -o "results/hqslq.csv" -a --format=$lineFormat timeout $timeLength solvers/HQSnp $pathToFile > "results/hqslqOutput/$fileName"
done