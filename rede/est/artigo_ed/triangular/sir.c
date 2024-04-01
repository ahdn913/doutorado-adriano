// gcc violence_rm.c -O3 -Wall -lgsl -lgslcblas -lm -o violence_rm
// ./violence_rm
// for I in $(seq 10000); do echo $I; ./violence $I ;done
// grep "^1" dat_b075/10000_30/percolacao.dat | wc -l
// 1 = Moderados, 2 = Radicais

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 100   // L = Ni = Nj
#define Ni 100  // tamanho da rede em i
#define Nj 100  // tamanho da rede em j
#define bet 0.7 // S + I -> I + I
#define gam 0.3 // I -> R
#define inf 0.005 // densidade inicial de infectados
#define Ne 3    // número de classes de indivíduos

void op(int t, int *phi){
	int i, j;
	FILE *arq1;
	char nome4[100];

	sprintf(nome4, "dat_snap/%d/vio-%d.dat", L, t);
	arq1= fopen(nome4, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq1, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq1, "\n");
	}
    fclose(arq1);
}

int main(int argc, char **argv){
	int i, j, x, t, nx, vx; 
	int ativo, passivo, vizinho, pos;
	int *phi;
	char nome2[100];
	char nome3[100];
	char nome5[100];
	double acao, ic, p_del, p_gam, p_eps;
	FILE *arquivo;
	FILE *arq2;
	FILE *arq5;
		
	phi= (int *) calloc(L*L, sizeof(int));
	
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(i= 0; i< L; i++){
		for(j= 0; j< L; j++){
			ic= gsl_rng_uniform(w);
			if(ic > rad){
				phi[i*Nj+j]= 1;
			}
			else{
				phi[i*Nj+j]= 2;
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
	sprintf(nome5, "curve/%d/data.dat", L);
	arq5= fopen(nome5, "w");
	while(t <= 100){
		fprintf(arq5, "%d ", t);
		nx= 1;
		while(nx != Ne+1){
			vx= 0;
			for(i= 0; i< Ni; i++){
				for(j= 0; j< Nj; j++){
					if(phi[i*Nj+j] == nx){
						vx++;
					}
				}
			}
			fprintf(arq5, "%d ", vx);
			if(nx == Ne){
				fprintf(arq5, "\n");
			}
			nx++;
		}
		for(x= 0; x< Ni*Nj; x++){
			i= gsl_rng_uniform(w)*Ni;
			j= gsl_rng_uniform(w)*Nj;
			ativo= i*Nj+j;
			if(phi[ativo] != 0){
				if(i%2 == 0){
					vizinho= gsl_rng_uniform(w)*6;
					switch(vizinho){
						case 0:
							passivo= i*Nj+j+1;
							pos= 1; 	// right
						break;
						case 1:
							passivo= i*Nj+j-L;
							pos= 2; 	// top right
						break;
						case 2:
							passivo= i*Nj+j-L-1;
							pos= 3; 	// top left
						break;
						case 3:
							passivo= i*Nj+j-1;
							pos= 4; 	// left
						break;
						case 4:
							passivo= i*Nj+j+L-1;
							pos= 5; 	// bottom left
						break;
						default:
							passivo= i*Nj+j+L;
							pos= 6; 	// bottom right
						break;
					}
					if(i != 0 && j == 0){ // borda esquerda
						if(pos == 3){
							passivo= i*Nj+j-1;
						}
						else if(pos == 4){
							passivo= i*Nj+j+L-1;
						}
						else if(pos == 5){
							passivo= i*Nj+j+(2*L)-1;
						}
					}
					else if(i != 0 && j == (L-1)){ // borda direita
						if(pos == 1){
							passivo= i*Nj+j-L+1;
						}
					}
					else if(i == 0 && j != 0 && j != (L-1)){ // borda superior
						if(pos == 2){
							passivo= i*Nj+j+(L*L)-L;
						}
						else if(pos == 3){
							passivo= i*Nj+j+(L*L)-L-1;
						}
					}
					else if(i == 0 && j == 0){ // borda canto superior esquerdo
						if(pos == 2){
							passivo= (L*L)-L;
						}
						else if(pos == 3){
							passivo= (L*L)-1;
						}
						else if(pos == 4){
							passivo= L-1;
						}
						else if(pos == 5){
							passivo= (2*L)-1;
						}
					}
					else if(i == 0 && j == (L-1)){ // borda canto superior direito
						if(pos == 1){
							passivo= 0;
						}
						else if(pos == 2){
							passivo= (L*L)-1;
						}
						else if(pos == 3){
							passivo= (L*L)-2;
						}
					}
				}
				else{
					vizinho= gsl_rng_uniform(w)*6;
					switch(vizinho){
						case 0:
							passivo= i*Nj+j+1;
							pos= 1;
						break;
						case 1:
							passivo= i*Nj+j-L+1;
							pos= 2;
						break;
						case 2:
							passivo= i*Nj+j-L;
							pos= 3;
						break;
						case 3:
							passivo= i*Nj+j-1;
							pos= 4;
						break;
						case 4:
							passivo= i*Nj+j+L;
							pos= 5;
						break;
						default:
							passivo= i*Nj+j+L+1;
							pos= 6;
						break;
					}
					if(j == 0 && i != (L-1)){ // borda esquerda
						if(pos == 4){
							passivo= i*Nj+j+L-1;
						}
					}
					else if(j == (L-1) && i != (L-1)){ // borda direita
						if(pos == 1){
							passivo= i*Nj+j-L+1;
						}
						else if(pos == 2){
							passivo= i*Nj+j-(2*L)+1;
						}
						else if(pos == 6){
							passivo= i*Nj+j+1;
						}
					}
					else if(i == (L-1) && j != 0 && j != (L-1)){ // borda inferior
						if(pos == 5){
							passivo= (i*Nj+j)%L;
						}
						if(pos == 6){
							passivo= ((i*Nj+j)%L)+1;
						}
					}
					else if(i == (L-1) && j == 0){ // borda canto inferior esquerdo
						if(pos == 4){
							passivo= (L*L)-1;
						}
						else if(pos == 5){
							passivo= 0;
						}
						else if(pos == 6){
							passivo= 1;
						}
					}
					else if(i == (L-1) && j == (L-1)){ // borda canto inferior direito
						if(pos == 1){
							passivo= (L*L)-L;
						}
						else if(pos == 2){
							passivo= (L*L)-2*L;
						}
						else if(pos == 5){
							passivo= L-1;
						}
						else if(pos == 6){
							passivo= 0;
						}
					}
				}
				
				acao= gsl_rng_uniform(w);

				if(acao < p_gam){ // ação de contágio de radical
					if(phi[ativo] == 2 && phi[passivo] == 1){
						phi[ativo]= phi[passivo]= 2;
					}
				}
				else if(acao < (p_gam + p_del)){ // ação de contágio de midias sociais
					if(phi[ativo] == 1){
						phi[ativo]= 2;
					}
				}
				else{
					if(phi[ativo] == 2){ // ação da coerção policial/judicial
						phi[ativo]= 1;
					}
				}	
			}
		}

		t++;
		op(t, phi);
		if(t%10 == 0){
			printf("%d \n", t);
		}

		sprintf(nome2, "dat_class/%d/%d.dat", L, t); // este bloco serve para criar os arquivos e, se já existirem, os deletarem.
		arquivo= fopen(nome2, "w");
		fclose(arquivo);

		sprintf(nome2, "dat_class/%d/%d.dat", L, t); // este bloco serve para escrever a quantidade de indivíduos nas espécies nos novos arquivos em branco
		arquivo= fopen(nome2, "a");
		fprintf(arquivo, "%d ", t);
		nx= 1;
		while(nx != Ne+1){
			vx= 1;
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
	fclose(arq5);

	gsl_rng_free(w);
	free(phi);
	return 0;
}
