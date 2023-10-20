// gcc -O3 -Wall percolacao.c -lgsl -lgslcblas -lm -o percolacao
// argv[1] = seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 1000 // L = Ni = Nj
#define Ni 1000   // tamanho da rede em i
#define Nj 1000   // tamanho da rede em j

int phi[L*L];

void op(int t){
	int i, j;
	FILE *arq;
	char nome[100];

	sprintf(nome, "dat/sir-%d.dat", t);
	arq= fopen(nome, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq, "\n");
	}
    fclose(arq);
}

void recur(int j, int i, int v) {
	if(i >= 0 && j >= 0 && i < L && j < L && phi[j*L+i] == 2){ 
		phi[j*L+i] = v;
		recur((j+1), i, v);
		recur((j-1), i, v);
		recur(j, (i+1), v);
		recur(j, (i-1), v);
	}
}

int count_clusters(void){
	int i, j, cls;
	cls= L*L; 
	for(j= 0; j< L; j++){
		for(i= 0; i< L; i++){
			if(phi[j*L+i] != 2) 
				continue ;
			recur(j, i, ++cls);
		}
	}
	return cls;
}

int main(int argc, char **argv){
	int i, j, x; 
	int pr_int;
	int max;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, lambda;
	FILE *file;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(lambda= 0.5854; lambda< 0.6010; lambda+= 0.0012){
		for(i= 0; i< Ni; i++){
			for(j= 0; j< Nj; j++){
				acao= gsl_rng_uniform(w);
				if(acao < lambda){
					phi[i*Nj+j]= 2;
				}
				else{
					phi[i*Nj+j]= 1;
          	  	}
			}
		}
		op(0);
		pr_int= lambda*10000;
		op(1);
	//				Percolação
		key= 0;
		max= count_clusters();
		for(i= L*L; i<= max; i++){ 
			a_ih= a_fh= 0;
			for(j= 0; j< L; j++){
				if(phi[j] == i){ 
					a_ih++; 
				}
				if(phi[(L-1)*L+j] == i){ 
					a_fh++; 
				}
			}

			if(a_ih > 0 && a_fh > 0){
				key= 1;
			}
		}
		
		sprintf(name, "dat/L%d/percolacao_%d.dat", L, pr_int);
		if(!(file= fopen(name, "a"))){
			printf("cannot open file teste\n");
			exit(0);
		}
		if(key == 1){
			fprintf(file, "1\n");
			//printf("percolou\n");
		}
		else{
			fprintf(file, "0\n");
			//printf("não percolou\n");
		}
		fclose(file);
	}
	gsl_rng_free(w);
	return 0;
}
