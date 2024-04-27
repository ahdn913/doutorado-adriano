// gcc sir.c -O3 -Wall -lgsl -lgslcblas -lm -o sir
// ./sir
// for I in $(seq 1 100); do echo $I; ./sir $I ;done
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado, 4 = Vacinado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 64 // L = Ni = Nj
#define Ni 64   // tamanho da rede em i
#define Nj 64   // tamanho da rede em j
#define NG 10000
//#define T 5.0
#define minT 0.5
#define difT 0.1
#define transiente 1000
#define norm 1.0/(NG*L*L)


int phi[L*L];

void op(int t){
	int i, j;
	FILE *arq;
	char nome[100];

	sprintf(nome, "dat_snap/%d/data-%d.dat", L, t);
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
	int i, j;
	double act;
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
			act= gsl_rng_uniform(w);
            if(act >= 0.5){
                phi[i*Nj+j]= 1;
            }
            else{
                phi[i*Nj+j]= -1;
            }
		}
	}
}

int main(int argc, char **argv){
	int i, j, x, a, b, A, B; 
	int ativo, passivo, vizinho, t=0;
	int baix, cima, dire, esqu, T10;
	char name[100];
	double acao, E, dE, M, Mabs, E_avg, etot, etotsq, Msq, Msq_avg, M_avg, mtot, mtotsq;
	double Mabs_avg, Mq_avg, mabstot, mqtot, EE=0, MM=0, Esq_avg, T=5.0;
	FILE *arq;
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus); 

	ic(phi);
	op(0);
	while(T>= minT){
		dE= 0;
		for(a= 0; a< transiente; a++){
			for(b= 0; b< L*L; b++){
				i= gsl_rng_uniform(w)*Ni;
				j= gsl_rng_uniform(w)*Nj;
				ativo= i*Nj+j;
				baix= ((i+1)%Ni)*Nj+j; // baixo
				cima= ((i-1+Ni)%Ni)*Nj+j; // cima
				dire= i*Nj+(j+1)%Nj; // direita
				esqu= i*Nj+(j-1+Nj)%Nj; // esquerda
				acao= gsl_rng_uniform(w);
				E= -1*phi[ativo]*(phi[esqu] + phi[dire] + phi[cima] + phi[baix]);
				dE= -2*E;
				if(dE < 0){
					// flipa devido à baixa energia
					phi[ativo]= -phi[ativo];
				}
				else if(acao < exp(-dE/T)){
					// flipa devido ao banho térmico
					phi[ativo]= -phi[ativo];
				}
			}
		}
	
		M= 0;
		for(i= 0; i < L; i++){
			for(j= 0; j < L; j++){
				ativo= i*Nj+j;
				M+= phi[ativo]; // M = magnetização total da rede
			}
		}
		Mabs= abs(M);

		E= 0;
		for(i= 0; i < L; i++){
			for(j= 0; j < L; j++){
				ativo= i*Nj+j;
				baix= ((i+1)%Ni)*Nj+j; // baixo
				cima= ((i-1+Ni)%Ni)*Nj+j; // cima
				dire= i*Nj+(j+1)%Nj; // direita
				esqu= i*Nj+(j-1+Nj)%Nj; // esquerda
				E+= -1*phi[ativo]*(phi[esqu] + phi[dire] + phi[cima] + phi[baix]);
			}
		}

		etot= 0;
		etotsq= 0;
		mtot= 0;
		mtotsq= 0;
		mabstot= 0;
		mqtot= 0;

		for(a= 0; a < NG; a++){
			for(b= 0; b < L; b++){
				i= gsl_rng_uniform(w)*Ni;
				j= gsl_rng_uniform(w)*Nj;
				ativo= i*Nj+j;
				baix= ((i+1)%Ni)*Nj+j; // baixo
				cima= ((i-1+Ni)%Ni)*Nj+j; // cima
				dire= i*Nj+(j+1)%Nj; // direita
				esqu= i*Nj+(j-1+Nj)%Nj; // esquerda
				acao= gsl_rng_uniform(w);
				E= -1*phi[ativo]*(phi[esqu] + phi[dire] + phi[cima] + phi[baix]);
				dE= -2*E;
				if(dE < 0){
					// flipa devido à baixa energia
					phi[ativo]= -phi[ativo];
					EE+= 2*dE;
					MM+= 2*phi[ativo];
					Mabs+= abs(phi[ativo]);
				}
				else if(acao < exp(-dE/T)){
					// flipa devido ao banho térmico
					phi[ativo]= -phi[ativo];
					EE+= 2*dE;
					MM+= 2*phi[ativo];
					Mabs+= abs(phi[ativo]);
				}
			}
			etot+= EE/2.0;
			etotsq+= (EE/2.0)*(EE/2.0);
			mtot+= MM;
			mtotsq+= MM*MM;
			mqtot+= MM*MM*MM*MM;
			mabstot+= (sqrt(M*M));
		}

		E_avg= etot*norm;
		Esq_avg= etotsq*norm;
		M_avg= mtot*norm;
		Msq_avg= mtotsq*norm;
		Mabs_avg= mabstot*norm;
		Mq_avg= mqtot*norm; 
	
		T10= T*10;
		sprintf(name, "dat_snap/%d/data-%d.dat", L, T10);
		arq= fopen(name, "w");
		fprintf(arq, "<M> <|M|> <M^2> X X' <E> <E^2> C U_L\n");
		fprintf(arq, "%e %e %e %e %e %e %e %e %e\n", M_avg, Mabs_avg, Msq_avg, (Msq_avg-(M_avg*M_avg*L*L)/T), (Msq_avg-(Mabs_avg*Mabs*L*L)/T), E_avg, Esq_avg, Esq_avg-(E_avg*E_avg*L*L)/(T*T), (1-((Mq_avg)/(3*Msq_avg))));
		fclose(arq);
		T+=-difT;
	}

	gsl_rng_free(w);
	return 0;
}
