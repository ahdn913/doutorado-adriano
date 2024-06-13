// gera a posicao dos individuos no espaco
// gcc x_y.c -O3 -Wall -lgsl -lgslcblas -lm -o x_y
// ./x_y numero_de_individuos seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

int main(int argc, char **argv){
    int i, N;
    double x, y;
    FILE *arq;
    char name[100];

    N= atoi(argv[1]);

    gsl_rng_default_seed = atoi(argv[2]);
	gsl_rng *w = gsl_rng_alloc(gsl_rng_taus);

    sprintf(name, "xy/xy-%d-%d.dat", atoi(argv[1]), atoi(argv[2]));
	arq = fopen(name, "w");

    for(i= 0; i< N; i++){
        x= gsl_rng_uniform(w);
        y= gsl_rng_uniform(w);
        fprintf(arq, "%e %e\n", x, y);
    }
    fclose(arq);
    
    gsl_rng_free(w);
    return 0;
}