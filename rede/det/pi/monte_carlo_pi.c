// Modelo para calcular pi
// compilação: gcc monte_carlo_pi.c -O3 -Wall -lgsl -lgslcblas -lm -o monte_carlo_pi
// for I in $(seq 10000); do echo $I; ./monte_carlo_pi $I ;done
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define N 10000  

struct intble {
	double x; 
	double y;
};

int main(int argc, char **argv){
	int i, cont_d=0, seed;
	double d, pi;
	struct intble *p;
	FILE *arq;
	char clus[100];

	seed= atoi(argv[1]); 

	p= (struct intble *) calloc(N*N, sizeof(struct intble)); // Calloc aloca uma matriz na memória (número, tamanho)

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); // inicialização do gerador de número aleatório (algorítmo Taus)
	for(i= 0; i< N; i++){					
		p[i].x= gsl_rng_uniform(w);
		p[i].y= gsl_rng_uniform(w);
	}

	for(i= 0; i< N; i++){
		d= sqrt(p[i].x*p[i].x + p[i].y*p[i].y);
		if(d <= 1.0){
			cont_d++;
		}
	}
	pi= 4.0*cont_d/N;
	sprintf(clus, "dat/pi-%d.dat", seed); 
	arq= fopen(clus, "a");
	fprintf(arq, "%e\n", pi); 
	fclose(arq);
	gsl_rng_free(w);
	free(p);
	return 0;
}