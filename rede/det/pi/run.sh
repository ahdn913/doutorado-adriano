#!/bin/bash

gcc monte_carlo_pi.c -O3 -Wall -lgsl -lgslcblas -lm -o monte_carlo_pi
gcc avg.c -O3 -Wall -lgsl -lgslcblas -lm -o avg
for J in $(seq 1000)
do
    echo $J
    ./monte_carlo_pi $J
done
./avg pi 1 1 1000
