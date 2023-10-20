// gcc sir_vax_b075_otim.c -O3 -Wall -lgsl -lgslcblas -lm -o sir_vax_b075_otim
// ./sir_vax_b075_otim
// for I in $(seq 10000); do echo $I; ./sir_vax_b075_otim $I ;done
// grep "^1" dat_b075/10000_30/percolacao.dat | wc -l
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado, 4 = Vacinado
// argv[1] = seed
//teste

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 30 // L = Ni = Nj
#define Ni 30   // tamanho da rede em i
#define Nj 30   // tamanho da rede em j
#define repeticoes 10000

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

void ic(int *phi){ 
	int i, j, cont_inf;

	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
            if(i == Ni/2 && j == Nj/2){
                phi[i*Nj+j]= 2;
				cont_inf= 1;
            }
            else{
                phi[i*Nj+j]= 1;
            }
		}
	}
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
	int ativo, passivo, vizinho, pr_int, cont_inf, cont_sus, cont_otim, geracao;
	int max;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, lambda, vax, epsilon;
	FILE *file;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(epsilon= 0.0222; epsilon< 0.0282; epsilon+= 0.0002){
		lambda=1.0;
		ic(phi);
		cont_sus= 1;
		op(0);
		pr_int= epsilon*1000000; 
		geracao= 0;
		while(cont_sus > 0){
			for(x= 0; x< Ni*Nj; x++){
				i= gsl_rng_uniform(w)*Ni;
				j= gsl_rng_uniform(w)*Nj;
				ativo= i*Nj+j;
				vizinho= gsl_rng_uniform(w)*4;
				switch(vizinho){
					case 0:
						passivo= ((i+1)%Ni)*Nj+j;
					break;
					case 1:
						passivo= ((i-1+Ni)%Ni)*Nj+j;
					break;
					case 2:
						passivo= i*Nj+(j+1)%Nj;
					break;
					default:
						passivo= i*Nj+(j-1+Nj)%Nj;
					break;
				}
				acao= gsl_rng_uniform(w);
				if(acao < lambda){
					//infecção
					if(phi[passivo] == 1 && phi[ativo] == 2){
						phi[passivo]= phi[ativo];
					}
				}
				/*else{
					//recuperação
					if(phi[ativo] == 2){
						phi[ativo]= 3;
					}
				}*/
				if(phi[ativo]==1){
					vax= gsl_rng_uniform(w);
					if(epsilon > vax){
						phi[ativo]=4;
					}
				}
			}
			cont_sus= 0;
			for(i= 0; i < Ni; i++){
				for(j= 0; j < Nj; j++){
					if(phi[i*Nj+j]==1){
						cont_sus++;
					}
				}
			}
			geracao++;
			cont_otim= 0;
			if(geracao % 200 == 0){
				for(i= 0; i < Ni; i++){
					for(j= 0; j < Nj; j++){
						if(phi[i*Nj+j] == 2){
							if(phi[((i+1)%Ni)*Nj+j] == 1){
								cont_otim++;
							}
							if(phi[((i-1+Ni)%Ni)*Nj+j] == 1){
								cont_otim++;
							}
							if(phi[i*Nj+(j+1)%Nj] == 1){
								cont_otim++;
							}
							if(phi[i*Nj+(j-1+Nj)%Nj] == 1){
								cont_otim++;
							}
						}
					}
				}
				if(cont_otim == 0){
					cont_sus= 0;
				}
			}
			// adicionar os 2 loop for para transformar todos os suscetíveis em vacinados se cont_otim == 0?
		}
		op(1);
	//				Percolação
		key= 0;
		max= count_clusters();
		for(i= L*L+1; i<= max; i++){ 
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
		
		sprintf(name, "dat_b075/%d_%d/percolacao_%d.dat", repeticoes, L, pr_int);
		if(!(file= fopen(name, "a"))){
			printf("cannot open file teste\n");
			exit(0);
		}
		if(key == 1){
			fprintf(file, "1\n");
			printf("percolou\n");
		}
		else{
			fprintf(file, "0\n");
			printf("não percolou\n");
		}
		fclose(file);
	}
	gsl_rng_free(w);
	return 0;
}
