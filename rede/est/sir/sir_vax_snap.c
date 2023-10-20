// gcc sir_vax_snap.c -O3 -Wall -lgsl -lgslcblas -lm -o sir_vax_snap
// ./sir_vax_snap
// for I in $(seq 100 199); do echo $I; ./sir_vax_snap $I ;done
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
#define repeticoes 10000

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
	int ativo, passivo, vizinho, pr_int, cont_inf, cont_sus, cont_s, cont_i, cont_r, cont_v;
	int max, t=0;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, lambda, vax, epsilon;
	FILE *file;
	FILE *arq;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(epsilon= 0.0020; epsilon< 0.0021; epsilon+= 0.0002){
		lambda=0.85;
		ic(phi);
		cont_sus= 1;
		cont_inf= 1;
		pr_int= epsilon*10000000; //*************** diminuir 1 zero
		op(0);
		while(cont_sus > 0){
			while(cont_inf > 0){
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
				sprintf(name, "plot_curve/data-%d.dat", L); // Retorna uma string de acordo com a string de formatação Format
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
					if(acao < lambda){
						//infecção
						if(phi[passivo] == 1 && phi[ativo] == 2){
							phi[passivo]= phi[ativo];
						}
					}
					else{
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
				if(cont_inf == 0){
					cont_sus= 0;
					t++;
					sprintf(name, "plot_curve/%d/data.dat", L); // Retorna uma string de acordo com a string de formatação Format
					arq= fopen(name, "a"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
					fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
					fclose(arq);
					t+=-1;
				}
				t++;
				op(t);
			}
		}

		
	//				Percolação
		/*
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
		
		sprintf(name, "dat/%d_%d/percolacao_%d.dat", repeticoes, L, pr_int);
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
		*/
	}
	gsl_rng_free(w);
	return 0;
}
