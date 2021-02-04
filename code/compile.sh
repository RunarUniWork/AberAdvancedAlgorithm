#!/bin/bash
for each in $(ls *.c)
do	
	echo "${each}"
	gcc -Wall $each -lglpk
	#./a.out inputs/medium
	./a.out -d inputs/small222
done
