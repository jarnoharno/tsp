#/bin/bash
for m in greedy mstp cw
do
	for f in input/*
	do
		./tsp $m $f > out/`basename $f`.${m}.out
	done
done
