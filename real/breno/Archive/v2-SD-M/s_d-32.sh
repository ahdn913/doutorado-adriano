#!/bin/bash
 
# created     : 2024/06/02
# last update : 2024/06/02
# author      : breno <bfoliveira@uem.br>
# notes       : 

if [[ -f dat/S_T_C-32-1.dat  ]]
then
	N=$(wc -l dat/S_T_C-32-1.dat | rev | cut -f2 -d" " | rev)
	S0=$(echo "$N/101 - 100" | bc)
	T0=$(echo "$N%101 + 100" | bc)
fi

for S in $(seq $S0 2 100)
do
	for T in $(seq $T0 2 200)
	do
		./s_d.out inp/n_b-250000-32-13.dat inp/x_y-250000-32-13.dat 250000 $S $T 1 32
	done
done
