#!/bin/bash

gcc -Wall puzzleSolver.c -lglpk

for each in $(ls inputs/in*.txt)
do	
	#echo "${each}"
	./a.out ${each}
	echo ""
	echo "================="
	echo ""
done
