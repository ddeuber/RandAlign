#!/bin/bash

KS="8 64 512 4096 32768"

for K in $KS; do 
    export PARALLEL_READS=${K}
    time(./app run 'data/genome.chr22.fa.txt' 'data/output_5xCov1.fq' 'data/output_5xCov2.fq'> /dev/null) 2>> timings.txt
done 