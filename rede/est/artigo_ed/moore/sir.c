// gcc sir.c -O3 -Wall -lgsl -lgslcblas -lm -o sir
// ./sir
// for I in $(seq 10000); do echo $I; ./sir $I ;done
// grep "^1" dat_b075/10000_30/percolacao.dat | wc -l
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado
// argv[1] = seed
//teste

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 100 // L = Ni = Nj
#define Ni 100   // tamanho da rede em i
#define Nj 100   // tamanho da rede em j

int phi[L*L];

void op(int t){
	int i, j;
	FILE *arq;
	char nome[100];

	sprintf(nome, "dat_snap/vio-%d.dat", t); // só é usado para fazer os snapshots
	arq= fopen(nome, "w");
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			fprintf(arq, "%d ", phi[i*Nj+j]);
		}
		fprintf(arq, "\n");
	}
    fclose(arq);
}

// void ic(int *phi){ // condição inicial do sistema
// 	int i, j;
// 	double teste;

// 	gsl_rng *a= gsl_rng_alloc(gsl_rng_taus);
// 	teste= gsl_rng_uniform(a);

// 	for(i= 0; i< Ni; i++){
// 		for(j= 0; j< Nj; j++){
//             if(teste < 0.5){
//                 phi[i*Nj+j]= 2;
//             }
//             else{
//                 phi[i*Nj+j]= 1;
//             }
// 		}
// 	}
// }

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
	int i, j, x, xis; 
	int ativo, passivo, vizinho, pr_int, cont_inf, cont_sus, cont_otim, geracao, pr, pm, gamma1, t1;
	int max;
	int a_ih, a_fh, a_iv, a_fv, key;
	char name[100];
	double acao, gamma, rec, epsilon, teste, Pr;
	FILE *file;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 
	for(gamma= 0.7265; gamma<= 0.7268; gamma+= 0.0001){
		epsilon= 1-gamma;
		//op(0);
		t1= 1;
		while(t1 != 1001){
			for(i= 0; i< Ni; i++){
				for(j= 0; j< Nj; j++){
					teste= gsl_rng_uniform(w);
					if(teste < 0.5){
						phi[i*Nj+j]= 2;
					}
					else{
						phi[i*Nj+j]= 1;
					}
				}
			}
			if(t1 % 20 == 0){
				printf("Gamma atual: %e. simulações realizadas: %d\n", gamma, t1);
			}
			gamma1= gamma*10000;
			pr_int= gamma*10000; 
			sprintf(name, "dat_perco/%d/%d/perco_%d.dat", L, gamma1, t1);
			file= fopen(name, "w");
			fclose(file);
			geracao= 0;
			while(geracao <= 200){
				if(geracao >= 0){
					pr= pm= 0;
					Pr= 0.0;
					for(i= 0; i< Ni; i++){
						for(j= 0; j< Nj; j++){
							if(phi[i*Nj+j] == 2){
								Pr+= 1.0/(Ni*Nj);
							}
						}
					}
					//Pr= pr*1.0/(Ni*Nj*1.0);
					sprintf(name, "dat_perco/%d/%d/perco_%d.dat", L, gamma1, t1);
					file= fopen(name, "a");
					fprintf(file, "%d %e\n", geracao, Pr);
					fclose(file);
				}
				for(x= 0; x< Ni*Nj; x++){
					i= gsl_rng_uniform(w)*Ni;
					j= gsl_rng_uniform(w)*Nj;
					ativo= i*Nj+j;
					xis= i*Nj+j;
					vizinho= gsl_rng_uniform(w)*4;
					
					if(i!=0 && j!=0 && i!=(L-1) && j!=(L-1)){
						switch(vizinho){
							case 0:
								passivo= xis-L+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo=xis-L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo=xis+L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo=xis+L+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(j==0 && i!=0 && i!=(L-1)){ // esquerda
						switch(vizinho){
							case 0:
								passivo= xis-L+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo=xis-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo=xis+2*L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo=xis+L+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(j==(L-1) && i!=0 && i!=(L-1)){ // direita
						switch(vizinho){
							case 0:
								passivo= xis-2*L-1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo=xis-L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo=xis+L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo= xis+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==(L-1) && j!=0 && j!=(L-1)){ // baixo
						switch(vizinho){
							case 0:
								passivo= xis-L+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo=xis-L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo= xis-L*(L-1)-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo= xis-L*(L-1)+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==0 && j!=0 && j!=(L-1)){ // cima
						switch(vizinho){
							case 0:
								passivo= xis+L*(L-1)+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo= xis+L*(L-1)-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo=xis+L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo=xis+L+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==0 && j==0){ // superior esquerdo
						switch(vizinho){
							case 0:
								passivo= xis+L*(L-1)+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo= xis+L*L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo= xis+2*L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo=xis+L+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==0 && j==(L-1)){ // superior direito
						switch(vizinho){
							case 0:
								passivo= L*L-L;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo= L*L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo=xis+L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo= xis+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==(L-1) && j==(L-1)){
						switch(vizinho){
							case 0:
								passivo= xis-2*L+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo=xis-L-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo= L-2;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo= xis+1;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					else if(i==(L-1) && j==0){
						switch(vizinho){
							case 0:
								passivo= xis-L+1;
							break;
							case 1:
								passivo= ((i-1+Ni)%Ni)*Nj+j;
							break;
							case 2:
								passivo= xis-1;
							break;
							case 3:
								passivo= i*Nj+(j-1+Nj)%Nj;
							break;
							case 4:
								passivo= L-1;
							break;
							case 5:
								passivo= ((i+1)%Ni)*Nj+j;
							break;
							case 6:
								passivo=xis+1-(L*L)+L;
							break;
							case 7:
								passivo= i*Nj+(j+1)%Nj;
							break;
						}
					}
					


				// switch(vizinho){
				// case 0:
				// 	passivo= ((i+1)%Ni)*Nj+j;
				// break;
				// case 1:
				// 	passivo= ((i-1+Ni)%Ni)*Nj+j;
				// break;
				// case 2:
				// 	passivo= i*Nj+(j+1)%Nj;
				// break;
				// default:
				// 	passivo= i*Nj+(j-1+Nj)%Nj;
				// break;
				// }	

					acao= gsl_rng_uniform(w);
					if(acao < gamma){
						//infecção
						if(phi[passivo] == 1 && phi[ativo] == 2){
							phi[passivo]= 2;
						}
					}
					else{
						phi[ativo]=1;
					}
					
				}
				geracao++;
			}
			t1++;
		}
		//op(1);
	/*//				Percolação
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
		
		sprintf(name, "dat_perco/%d/perco_%d.dat", L, pr_int);
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
