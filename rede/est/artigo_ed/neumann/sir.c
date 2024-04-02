// gcc sir.c -O3 -Wall -lgsl -lgslcblas -lm -o sir
// ./sir
// for I in $(seq 1 100); do echo $I; ./sir $I ;done
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado, 4 = Vacinado
// argv[1] = seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 200 // L = Ni = Nj
#define Ni 200   // tamanho da rede em i
#define Nj 200   // tamanho da rede em j

int phi[L*L];

void op(int t){
	int i, j;
	FILE *arq;
	char nome[100];

	sprintf(nome, "dat_snap/%d/sir-%d.dat", L, t);
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
	int ativo, passivo, vizinho, cont_inf, cont_sus, cont_s, cont_i, cont_r, cont_v, trigger;
	int max, t=0;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, lambda, gamma, vax, epsilon, acao2;
	FILE *file;
	FILE *arq;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	sprintf(name, "curve/data.dat"); // Retorna uma string de acordo com a string de formatação Format
	arq= fopen(name, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
	fclose(arq);

	trigger= 0;
	lambda= 0.85;
	gamma= 0.15;
	ic(phi);
	cont_sus= 1;
	cont_inf= 1;
	op(0);
	while(cont_sus > 0){
		while(cont_inf > 0){
			while(t < 500){
				cont_s = cont_i = cont_r = cont_v = 0;
				for(i= 0; i< L*L; i++){
					if(phi[i] == 1){
						cont_s++;
					}
					else if(phi[i] == 2){
						cont_i++;
					}
					else if(phi[i] == 3){
						cont_r++;
					}
					else if(phi[i] == 4){
						cont_v++;
					}
				}
				sprintf(name, "curve/data.dat"); // Retorna uma string de acordo com a string de formatação Format
				arq= fopen(name, "a"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
				fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
				fclose(arq);

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
					acao2= gsl_rng_uniform(w);
					if(acao < lambda){
						//infecção
						if(phi[passivo] == 1 && phi[ativo] == 2){
							phi[passivo]= phi[ativo];
						}
					}
					if(acao2 < gamma){
						//recuperação
						if(phi[ativo] == 2){
							phi[ativo]= 3;
						}
					}
					/*if(phi[ativo]==1){
						vax= gsl_rng_uniform(w);
						if(epsilon > vax){
							phi[ativo]=4;
						}
					}*/
				}
				cont_sus= 0;
				for(i= 0; i < Ni; i++){
					for(j= 0; j < Nj; j++){
						if(phi[i*Nj+j]==1){
							cont_sus++;
						}
					}
				}
				cont_inf= 0;
				for(i= 0; i < Ni; i++){
					for(j= 0; j < Nj; j++){
						if(phi[i*Nj+j]==2){
							cont_inf++;
						}
					}
				}
				if(cont_inf == 0 && trigger == 0){
					printf("O último infectado foi recuperado em t = %d. \n", t-1);
					trigger= 1;
				}
				if(cont_inf == 0){
					if(t>30){
						cont_sus= 0;
						t++;
						sprintf(name, "curve/data.dat"); // Retorna uma string de acordo com a string de formatação Format
						arq= fopen(name, "a"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
						t+=-1;
					}
					else{
						ic(phi);
						t= 0;
						sprintf(name, "curve/data.dat"); // Retorna uma string de acordo com a string de formatação Format
						arq= fopen(name, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
					}
				}
				t++;
				op(t);
			}
		}
	}
	gsl_rng_free(w);
	return 0;
}
