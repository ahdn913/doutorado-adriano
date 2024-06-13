#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NG 10000
#define K 0.1
int main(int argc, char **argv){
	if(argc != 8){
		printf("%s n_b.dat x_y.dat N S T seed M\n", argv[0]);
		exit(1);
	}

	int i, j, k, t, ac, nb, n_b_max;
	int *p, *n_b;
	double C, S, T, P_ac, P_nb;
	char name[100];
	FILE *file;

	int N = atoi(argv[3]);
	S = atoi(argv[4])/100.0;
	T = atoi(argv[5])/100.0;
	const double payoff[4]= {\
		0.0, T,   \
		S,   1.0, \
	};
	if((p = (int *) calloc(N, sizeof(int))) == NULL){
		printf("cannot allocate p\n");
		exit(1);
	}
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
			if((fscanf(file, "%d", &ac) != 1)){
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

	gsl_rng_default_seed = atoi(argv[6]);
	gsl_rng *w = gsl_rng_alloc(gsl_rng_taus);
	for(i = 0; i < N; i++){
		if(gsl_rng_uniform(w) < 0.5){
			p[i] = 1; // cooperator
		}else{
			p[i] = 0; // desertor
		}
	}

	C = 0.0;
	for(t = 1; t <= NG; t++){
		for(k = 0; k < N; k++){
			ac = gsl_rng_uniform(w)*N;
			P_ac = 0.0;
			for(i = 1; i <= n_b[ac*n_b_max]; i++){
				P_ac+= payoff[p[ac]*2+p[n_b[ac*n_b_max+i]]];
			}
			nb = gsl_rng_uniform(w)*(n_b[ac*n_b_max])+1;
			nb = n_b[ac*n_b_max+nb];
			P_nb = 0.0;
			for(i = 1; i <= n_b[nb*n_b_max]; i++){
				P_nb+= payoff[p[nb]*2+p[n_b[nb*n_b_max+i]]];
			}
			if(gsl_rng_uniform(w) < 1.0/(1.0+exp((P_ac-P_nb)/K))){
				p[ac] = p[nb];
			}
		}
		if(t > 5000){
			for(i = 0; i < N; i++){
				if(p[i] == 1){
					C+= 1.0;
				}
			}
//			sprintf(name, "dat/d_C-%d-%d.dat", atoi(argv[7]), atoi(argv[6]));
//			file = fopen(name, "a");
//			fprintf(file, "%e\n", C/(N*(t-5000)));
//			fclose(file);
		}
	}
	free(n_b);
	sprintf(name, "dat/S_T_C-%d-%d.dat", atoi(argv[7]), atoi(argv[6]));
	file = fopen(name, "a");
	fprintf(file, "%e %e %e\n", S, T, C/(N*5000));
	fclose(file);

//	double *x, *y;
//	if(!(file = fopen(argv[2], "r"))){
//		printf("cannot open file %s\n", argv[2]);
//		exit(1);
//	}
//	if((x = (double *) calloc(N, sizeof(double))) == NULL){
//		printf("cannot allocate x\n");
//		exit(1);
//	}
//	if((y = (double *) calloc(N, sizeof(double))) == NULL){
//		printf("cannot allocate y\n");
//		exit(1);
//	}
//	for(i = 0; i < N; i++){
//		if((fscanf(file, "%lf %lf", &x[i], &y[i]) != 2)){
//			printf("cannot read file %s\n", argv[2]);
//			exit(1);
//		}
//	}
//	fclose(file);
//	sprintf(name, "dat/p-%d-%d.dat", atoi(argv[6]), atoi(argv[5]));
//	file = fopen(name,  "w");
//	for(i = 0; i< N; i++){
//		fprintf(file, "%d %.4e %.4e\n", p[i], x[i], y[i]);
//	}
//	fclose(file);
//	free(x);
//	free(y);

	free(p);
	gsl_rng_free(w);
	return 0;
}
