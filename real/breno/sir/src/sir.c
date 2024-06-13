#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NG 2000 // gerações
//#define NS 3 número de especies
#define pm 0.5 // mobilidade
/*#define PROBABILITY_OF_PREDATION \
const double pp[NS*NS]= {\
	0.0, 1.0, 0.0, \
	0.0, 0.0, 1.0, \
	1.0, 0.0, 0.0, \
};*/

int main(int argc, char **argv){
	if(argc != 6){
		printf("%s n_b.dat x_y.dat N M seed\n", argv[0]); // o n_b e o arquivo que contem o numero de vizinhos do indivíduo por linha e suas identidades
		exit(1);										  // x_y a posicao de cada um
	}

	int i, j, k, t, ac, gd, nb, tp, n_b_max;
	int *p, *d, *n_b, *prey;
	PROBABILITY_OF_PREDATION
	double at;
	double *x, *y;
	char name[100];
	FILE *file;

	int N = atoi(argv[3]);

	if(!(file = fopen(argv[1], "r"))){
		printf("cannot open file %s\n", argv[1]);
		exit(1);
	}
	n_b_max = 0;
	for(i = 0; i < N; i++){
		if((fscanf(file, "%d", &nb) != 1)){
			printf("1 - cannot read file %s\n", argv[1]);
			exit(1);
		}
		if(n_b_max < nb){
			n_b_max = nb;
		}
		for(j = 0; j < nb; j++){
			if((fscanf(file, "%d", &t) != 1)){
				printf("2 - cannot read file %s\n", argv[1]);
				exit(1);
			}
		}
	}
	fclose(file);
	n_b_max++;
	if((n_b = (int *) calloc(N*n_b_max, sizeof(int))) == NULL){
		printf("cannot allocate n_b\n");
		exit(1);
	}
	if(!(file = fopen(argv[1], "r"))){
		printf("cannot open file %s\n", argv[1]);
		exit(1);
	}
	for(i = 0; i < N; i++){
		if((fscanf(file, "%d", &n_b[i*n_b_max]) != 1)){
			printf("3 - cannot read file %s\n", argv[1]);
			exit(1);
		}
		k = n_b[i*n_b_max];
		for(j = 1; j < k+1; j++){
			if((fscanf(file, "%d", &n_b[i*n_b_max+j]) != 1)){
				printf("4 - cannot read file %s\n", argv[1]);
				exit(1);
			}
		}
	}
	fclose(file);

	if((prey = (int *) calloc(n_b_max, sizeof(int))) == NULL){
		printf("cannot allocate prey\n");
		exit(1);
	}

	gsl_rng_default_seed = atoi(argv[5]);
	gsl_rng *w = gsl_rng_alloc(gsl_rng_taus);

	if((p = (int *) calloc(N, sizeof(int))) == NULL){
		printf("cannot allocate p\n");
		exit(1);
	}

	if((d = (int *) calloc(NS, sizeof(int))) == NULL){
		printf("cannot allocate d\n");
		exit(1);
	}

	for(i = 0; i < N; i++){
		p[i]= gsl_rng_uniform(w)*NS+1;
	}

	if(!(file = fopen(argv[2], "r"))){
		printf("cannot open file %s\n", argv[2]);
		exit(1);
	}
	if((x = (double *) calloc(N, sizeof(double))) == NULL){
		printf("cannot allocate x\n");
		exit(1);
	}
	if((y = (double *) calloc(N, sizeof(double))) == NULL){
		printf("cannot allocate y\n");
		exit(1);
	}
	for(i = 0; i < N; i++){
		if((fscanf(file, "%lf %lf", &x[i], &y[i]) != 2)){
			printf("cannot read file %s\n", argv[2]);
			exit(1);
		}
	}
	fclose(file);

	for(t = 1; t <= NG; t++){
		gd = 0;
		while(gd < N){
			at = gsl_rng_uniform(w);
			ac = gsl_rng_uniform(w)*N;
			if(at < pm){
				nb = gsl_rng_uniform(w)*(n_b[ac*n_b_max])+1;
				nb = n_b[ac*n_b_max+nb];
				tp= p[nb];
				p[nb]= p[ac];
				p[ac]= tp;
				gd++;
			}else{
				j = 0;
				for(i = 1; i < n_b[ac*n_b_max]+1; i++){
					if(pp[(p[ac]-1)*NS+p[n_b[ac*n_b_max+i]]-1] == 1){
						prey[j] = n_b[ac*n_b_max+i];
						j++;
					}
					if(j != 0){
						nb = gsl_rng_uniform(w)*j;
						p[prey[nb]] = p[ac];
						gd++;
					}
				}
			}
		}
		if(t > 1000){
			for(i = 0; i < NS; i++){
				d[i] = 0;
			}
			for(i = 0; i < N; i++){
				d[p[i]-1]++;
			}
			sprintf(name, "dat-%d/d-%d.dat", atoi(argv[4]), atoi(argv[5]));
			file = fopen(name, "a");
			for(i = 0; i < NS; i++){
				fprintf(file, "%e ", (double)d[i]/N);
			}
			fprintf(file, "\n");
			fclose(file);
		}
	}
//	sprintf(name, "dat-%d/p-%d-%d.dat", atoi(argv[4]), 1, atoi(argv[5]));
//	file = fopen(name,  "w");
//	for(i = 0; i< N; i++){
//		fprintf(file, "%d %.4e %.4e\n", p[i], x[i], y[i]);
//	}
//	fclose(file);

	free(prey);
	free(n_b);
	free(x);
	free(y);
	free(p);
	gsl_rng_free(w);
	return 0;
}
