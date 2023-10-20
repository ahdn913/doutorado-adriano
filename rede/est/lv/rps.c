// gcc rps.c -O3 -Wall -lgsl -lgslcblas -lm -o rps
// ./rps
// for I in $(seq 10000); do echo $I; ./rps $I ;done
// grep "^1" dat_b075/10000_30/percolacao.dat | wc -l
// 0 = vazio, 1 = espécia A, 2 = Espécie B, 3 = Espécie C

// modifiquei os "free"

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 500 // L = Ni = Nj
#define Ni 500   // tamanho da rede em i
#define Nj 500   // tamanho da rede em j
#define mob 0.857 // mobilidade
#define pre 0.143 // predação
#define Ne 3    // número de espécies

//int phi[L*L];

void op(int t, int *phi){
	int i, j, geracao;
	FILE *arq1;
	char nome4[100];

	sprintf(nome4, "dat_snap/%d/%d/rps-%d.dat", Ne, L, t);
	arq1= fopen(nome4, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq1, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq1, "\n");
	}
    fclose(arq1);
	//free(arq1);
}

/*void ic(int *phi){ 
	int i, j, k1, l, t;
    double c1, c2, k2;
    char nome5[100];
    FILE *arq2;

    gsl_rng *p= gsl_rng_alloc(gsl_rng_taus); 
    for(i= 0; i< L; i++){
        for(j= 0; j< L; j++){
            l= 1;
            k2= gsl_rng_uniform(p)*Ne;
            k1= k2;
            while(l <= Ne){
                c1= 1.0*(l-1)/Ne;
                c2= 1.0*l/Ne;
                if(c1 <= k2/Ne){
                    if(k2/Ne < c2){
                        h= l;
                        phi[i*L+j]= h;
                    }
                }
                l++;
            }
        }
    }
    sprintf(nome3, "data.dat");
    arq2= fopen(nome3, "w");
    for(i= 0; i< L; i++){
		for(j= 0; j< L; j++){
			fprintf(arq2, "%d ", phi[i*L+j]);
		}
		fprintf(arq2, "\n");
	}
    fclose(arq2);

	gsl_rng_free(p);
	//free(arq3);
}*/

int main(int argc, char **argv){
	int i, j, x, t, g, seed, nx, vx, l, k1, h; 
	int ativo, passivo, vizinho;
	int *phi;
	char name[100];
	char nome[100];
	char nome2[100];
	char nome3[100];
	double acao, c1, c2, k2;
	FILE *arquivo;
	FILE *arq;
	FILE *arq2;
		
	//p= (struct intble *) calloc(N*N, sizeof(struct intble));
	phi= (int *) calloc(L*L, sizeof(int));
	
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(i= 0; i< L; i++){		// condição inicial
        for(j= 0; j< L; j++){
            l= 1;
            k2= gsl_rng_uniform(w)*Ne;
            k1= k2;
            while(l <= Ne){
                c1= 1.0*(l-1)/Ne;
                c2= 1.0*l/Ne;
                if(c1 <= k2/Ne){
                    if(k2/Ne < c2){
                        h= l;
                        phi[i*Nj+j]= h;
                    }
                }
                l++;
            }
        }
    }

    sprintf(nome3, "data.dat");		// este bloco serve para escrever o estado do sistema num .dat
    arq2= fopen(nome3, "w");
    for(i= 0; i< L; i++){
		for(j= 0; j< L; j++){
			fprintf(arq2, "%d ", phi[i*L+j]);
		}
		fprintf(arq2, "\n");
	}
    fclose(arq2);

	op(0, phi);
	t= 0;
	while(t <= 5000){
		for(x= 0; x< Ni*Nj; x++){
			i= gsl_rng_uniform(w)*Ni;
			j= gsl_rng_uniform(w)*Nj;
			ativo= i*Nj+j;
			if(phi[ativo] != 0){
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

				if(acao < mob){							 // mobilidade
					g= phi[ativo];
					phi[ativo]= phi[passivo];
					phi[passivo]= g;
				}

				else{ 
					if(phi[passivo] != phi[ativo]){         						// predação generalizada
						if(phi[ativo] != 1){
							if(phi[passivo] != phi[ativo]-1 && phi[passivo] != 0){
								phi[passivo]= phi[ativo];
							}
							else{
								x--;
							}
						}
						else{
							if(phi[passivo] != Ne && phi[passivo] != 0){
								phi[passivo]= phi[ativo];
							}
							else{
								x--;
							}
						}
					}
					else{
						x--;
					}
				}
			}
			else{
				x--;
			}
		}

		t++;
		op(t, phi);

		sprintf(nome2, "dat_class/%d/%d/%d.dat", Ne, L, t); // este bloco serve para criar os arquivos e, se já existirem, os deletarem.
		arquivo= fopen(nome2, "w");
		fclose(arquivo);

		sprintf(nome2, "dat_class/%d/%d/%d.dat", Ne, L, t); // este bloco serve para escrever a quantidade de indivíduos nas espécies nos novos arquivos em branco
		arquivo= fopen(nome2, "a");
		fprintf(arquivo, "%d ", t);
		nx= 0;
		while(nx != Ne+1){
			vx= 0;
			for(i= 0; i< Ni; i++){
				for(j= 0; j< Nj; j++){
					if(phi[i*Nj+j] == nx){
						vx++;
					}
				}
			}
			fprintf(arquivo, "%d ", vx);
			nx++;
		}
		fclose(arquivo);
	}

	gsl_rng_free(w);
	free(phi);
	//free(arquivo);
	//free(arq);
	return 0;
}
