#!/bin/sh

benchmarkFiles="../benchmarks/dqbf18/*"

timeLength="30s"

for pathToFile in $benchmarkFiles; do
    fileName=$(basename $pathToFile)
    echo $fileName
    lineFormat=$fileName",%e,%x"
    time -q -o "results/simplesolver.csv" -a --format=$lineFormat timeout $timeLength ./solver 0 $pathToFile
    time -q -o "results/treesolver.csv" -a --format=$lineFormat timeout $timeLength ./solver 1 $pathToFile
    time -q -o "results/hqs.csv" -a --format=$lineFormat timeout $timeLength ../hqs/hqs $pathToFile
done