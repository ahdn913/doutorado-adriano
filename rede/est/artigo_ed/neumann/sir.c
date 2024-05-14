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
#include <gsl/gsl_randist.h>

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
            if(i == (Ni+1)/6 && j == (Nj+1)/6){
                phi[i*Nj+j]= 2;
				cont_inf= 1;
            }
            else{
                phi[i*Nj+j]= 1;
            }
			if(i%40 == 0 || j%40 == 0 || i == Ni-1 || j == Nj-1){
				phi[i*Nj+j]= 0;
			}
		}
	}
}
/*
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
*/
int main(int argc, char **argv){
	int i, j, x; 
	int ativo, passivo, vizinho, cont_inf, cont_sus, cont_s, cont_i, cont_r, cont_v, trigger;
	int max, t=0, m;
	int a_ih, a_fh, a_iv, a_fv, key, u, ux, uy, trig1= 1, xx, yy, trigger2;
	char name[100];
	double acao, lambda, gamma, vax, epsilon, acao2, pp, alfa, qq;
	FILE *file;
	FILE *arq;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 
	//gsl_rng *r= gsl_rng_alloc(gsl_rng_taus);
	//gsl_rng *s= gsl_rng_alloc(gsl_rng_taus);

	sprintf(name, "curve/data-2g-%s.dat", argv[1]); // Retorna uma string de acordo com a string de formatação Format
	arq= fopen(name, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
	fclose(arq);

	trigger= 0;
	// lambda= 0.88/(1+trigger2);
	// gamma= 1 - lambda;
	// alfa= 0.01;
	ic(phi);
	cont_sus= 1;
	cont_inf= 1;
	op(0);
	//while(cont_sus > 0){
		//while(cont_inf > 0){
			while(t < 3000){
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
				sprintf(name, "curve/data-2g-%s.dat", argv[1]); // Retorna uma string de acordo com a string de formatação Format
				arq= fopen(name, "a"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
				fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
				fclose(arq);
				//u= gsl_ran_poisson(s, 2);
				//printf("O número aleatorio da distribuição é: %e\n", pp);
				//printf("Poisson: %d\n", u);

				for(x= 0; x< Ni*Nj; x++){
					i= gsl_rng_uniform(w)*Ni;
					j= gsl_rng_uniform(w)*Nj;
					ativo= i*Nj+j;
					xx= i;
					yy= j;
					//trig1= 1;
					
					//while(trig1 == 1){
					//ux= uy= 0;

					pp= gsl_ran_cauchy(w, 0.3);
					if(pp >= Ni || pp <= -Ni){
						while(pp >= Ni || pp <= -Ni){
							pp= gsl_ran_cauchy(w, 0.3);
						}
					}
					qq= gsl_ran_cauchy(w, 0.3);
					if(qq >= Nj || qq <= -Nj){
						while(qq >= Nj || qq <= -Nj){
							qq= gsl_ran_cauchy(w, 0.3);
						}
					}
/*
					pp= gsl_ran_poisson(w, 0.38);
					qq= gsl_ran_poisson(w, 0.38);
					acao= gsl_rng_uniform(w);
					if(acao < 0.5){
						pp= -1*pp;
					}
					acao= gsl_rng_uniform(w);
					if(acao < 0.5){
						qq= -1*qq;
					}
*/
					ux= pp;
					uy= qq;
					trigger2= 0;
					while(xx+ux >= Ni || xx+ux <= 0){
						pp= gsl_ran_cauchy(w, 0.3);
						ux= pp;
					}
					while(yy+uy >= Nj || yy+uy <= 0){
						qq= gsl_ran_cauchy(w, 0.3);
						uy= qq;
					}
					/*
					if(xx >= Ni){
						xx+= -Ni;
					}
					if(xx < 0){
						xx+= Ni;
					}
					yy+= uy;
					if(yy >= Nj){
						yy+= -Nj;
					}
					if(yy < 0){
						yy+= Nj;
					}*/
					
					if(ux > 0){
						for(m= xx; m< xx+ux; m++){
							i= m;
							passivo= i*Nj+j;
							if(phi[passivo] == 0){
								trigger2++;
							}
						}
					}
					else{
						for(m= xx; m> xx+ux; m--){
							i= m;
							passivo= i*Nj+j;
							if(phi[passivo] == 0){
								trigger2++;
							}
						}
					}
					if(uy > 0){
						for(m= yy; m< yy+uy; m++){
							j= m;
							passivo= i*Nj+j;
							if(phi[passivo] == 0){
								trigger2++;
							}
						}
					}
					else{
						for(m= yy; m> yy+uy; m--){
							j= m;
							passivo= i*Nj+j;
							if(phi[passivo] == 0){
								trigger2++;
							}
						}
					}
					passivo= (xx+ux)*Nj+(yy+uy);
					lambda= 0.8/((1+trigger2*4));
					gamma= 0.15;
					alfa= 0.02;
					//xx+= ux;
					
					//passivo= xx*Nj+yy;
					
					//printf("pp %e ux %d qq %e uy %d ativo %d passivo %d\n", pp, ux, qq, uy, ativo, passivo);
					

					/*vizinho= gsl_rng_uniform(w)*4;
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
					}*/

					acao= gsl_rng_uniform(w);
					acao2= gsl_rng_uniform(w);
					if(phi[ativo] == 2){
						if(acao < lambda){
						//infecção
							if(phi[passivo] == 1 && phi[ativo] == 2){
								phi[passivo]= phi[ativo];
							}
						}
						if(acao < gamma){
							//recuperação
							if(phi[ativo] == 2){
								phi[ativo]= 3;
							}
						}
					}
					if(phi[ativo] == 3){
						if(acao2 < alfa){
							phi[ativo]= 1;
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
				/*
				if(cont_inf == 0 && trigger == 0){
					printf("O último infectado foi recuperado em t = %d. \n", t-1);
					trigger= 1;
				}
				*/
				if(cont_inf == 0){
					if(t>30){
						cont_sus= 0;
						t++;
						sprintf(name, "curve/data-2g-%s.dat", argv[1]); // Retorna uma string de acordo com a string de formatação Format
						arq= fopen(name, "a"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
						if(trigger == 0){
							trigger= 1;
							printf("O tempo é %d\n", t);
						}
						t+=-1;
					}
					else{
						ic(phi);
						t= 0;
						sprintf(name, "curve/data-2g-%s.dat", argv[1]); // Retorna uma string de acordo com a string de formatação Format
						arq= fopen(name, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
						printf("resetou\n");
					}
				}
				t++;
				//op(t);
			}
			// tt++;
			// if(tt%5==0){
			// 	printf("tempo: %d\n", tt);
			// }
		//}
	//}
	gsl_rng_free(w);
	//gsl_rng_free(r);
	//gsl_rng_free(s);
	return 0;
}
