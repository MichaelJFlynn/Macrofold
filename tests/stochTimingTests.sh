#!/bin/sh
touch stochTimeTests.txt
rm stochTimeTests.txt
touch stochTimeTests.txt
echo "Filename\tLength\tNormal\tminus10\tminus8\tminus6\tminus4" >> stochTimeTests.txt
for((i = 100; i <= 4000;i=i+100))
do
    echo $i
    python random_seq.py $i > ./randseqs/${i}.seq;
    ../src/time_test ./randseqs/${i}.seq >> stochTimeTests.txt;
done
