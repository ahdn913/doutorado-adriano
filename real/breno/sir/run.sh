#!/bin/bash
 
# created     : 2024/05/27
# last update : 2024/05/27
# author      : breno <bfoliveira@uem.br>
# notes       : 


for I in $(seq 1001 10000)
do
	for J in 8 16 32 64 128 256 512
	do
 		./rps.out inp/n_b-250000-$J-13.dat inp/x_y-250000-$J-13.dat 250000 $J $I
	done
done 
