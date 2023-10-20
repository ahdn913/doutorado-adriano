// gcc sir_v.c -O3 -Wall -lgsl -lgslcblas -lm -o sir_v
// ./sir_v
// for I in $(seq 10000); do echo $I; ./sir_v $I ;done
// grep "^1" dat_b075/10000_30/percolacao.dat | wc -l
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado, 4 = Vacinado
// argv[1] = seed
//teste

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
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
struct stat st = {0};

void op(int t, double epsilon){
	int i, j, geracao, pr_int=epsilon*10000, l=L;
	FILE *arq1;
	char nome3[100];
	char *teste[100];
	//seed= atoi(argv[1]);

	sprintf(*teste, "dat_snap/%d/%d", l, pr_int); ///////////////////////////////////// o erro ta dentro desse void.
	if (stat(*teste, &st) == -1) { // se a pasta já existir, não faz nada. Caso contrário, é criada.
		mkdir(*teste, S_IRWXG); //  S_IRWXG é permissão para ler, escrever ou executar.
	}
	//mkdir(teste, 0777);
	sprintf(nome3, "dat_snap/%d/%d/siv-%d.dat", L, pr_int, t);
	arq1= fopen(nome3, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq1, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq1, "\n");
	}
    fclose(arq1);
}

/*void op_final(int seed){
	int i, j;
	FILE *arq;
	char nome[100];
	//seed= atoi(argv[1]);

	sprintf(nome, "dat/siv-%d.dat", seed);
	arq= fopen(nome, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq, "\n");
	}
    fclose(arq);
}*/

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

/*void recur(int j, int i, int v) {
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
}*/

int main(int argc, char **argv){
	int i, j, x, t, seed; 
	int ativo, passivo, vizinho, pr_int, cont_inf, cont_sus, cont_otim, geracao, class_s, class_i, class_v;
	int max;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	char nome[100];
	char nome2[100];
	double acao, lambda, vax, epsilon;
	FILE *file;
	FILE *arquivo;
	FILE *arq;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 
	seed= atoi(argv[1]); 

	for(epsilon= 0.0156; epsilon< 0.0157; epsilon+= 0.0200){
		lambda= 1.0;
		ic(phi);
		cont_sus= 1;
		op(0, epsilon);
		pr_int= epsilon*1000000; 
		geracao= 0;
		t= 0;
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
			t++;
			op(t, epsilon);
			cont_otim= 0;
			if(geracao % 1 == 0){                      // mudar aqui talvez?
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
					class_s = class_i = class_v = 0;
					for(i= 0; i < Ni; i++){
						for(j= 0; j < Nj; j++){
							if(phi[i*Nj+j] == 1){
								class_s++;
							}
							if(phi[i*Nj+j] == 2){
								class_i++;
							}
							if(phi[i*Nj+j] == 4){
								class_v++;
							}
						}
					}
					sprintf(nome2, "dat_class/%d/%d.dat", L, seed);
					arquivo= fopen(nome2, "w");
					fprintf(arquivo, "%d %d %d %d", geracao, class_s, class_i, class_v);
					fclose(arquivo);
				}
			}
		}
		sprintf(nome, "dat/siv-%d.dat", seed);
		arq= fopen(nome, "w");
		for(i= 0; i< Ni; i++){
			for(j= 0; j< Nj; j++){
				fprintf(arq, "%d ", phi[i*Nj+j]);
			}
			fprintf(arq, "\n");
		}
		fclose(arq);
	//				Percolação
		/*key= 0;
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
		fclose(file);*/
	}
	gsl_rng_free(w);
	return 0;
}
