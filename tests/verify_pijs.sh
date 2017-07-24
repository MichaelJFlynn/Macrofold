#!/bin/sh
for((i = 1; i <= 10;i=i+1))
do
    echo $i
    python random_seq.py 100 > ./randseqs/${i}.seq;
    ../src/hybrid-ss -D ./randseqs/${i}.seq >> dynamicThreshTest.txt;
done
