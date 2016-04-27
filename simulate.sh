#!/bin/bash
#PBS -l nodes=6:ppn=20
#cd $PBS_O_WORKDIR

for i in 0.10 0.14 0.20 0.30 0.45 0.67 1;
do
	for j in 0.04 0.05 0.06 0.07 0.08;
	do
		for k in 200 400 800;
		do
			python Penguin.py $k $j $i
		done
	done
done
