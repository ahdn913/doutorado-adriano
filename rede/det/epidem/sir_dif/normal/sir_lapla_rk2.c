// COMPILAÇÃO: gcc -O3 sir_lapla_rk2.c -Wall -lm -lgsl -lgslcblas -o sir_lapla_rk2
// EXECUÇÃO: ./sir_lapla_rk2
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define Nx 200  // tamanho da rede
#define Ny 200  // tamanho da rede
#define dt 0.1  // intervalo de tempo
#define D  1.0  // Parâmetro difusivo (mobilidade)
#define tf 300 // tempo final
#define NF 3000 // número de arquivos de saída
#define beta 0.85
#define gamma 0.15

void op(int n, double *S, double *I, double *R){ // printa a matriz de estado de cada geração
	int x, y;
	FILE *arq;
	char nome[100];
	sprintf(nome, "dat/sir-%d.dat", n);
	arq= fopen(nome, "w");

	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			fprintf(arq, "%e ", 1.0*S[x*Ny+y] + 2.0*I[x*Ny+y] + 3.0*R[x*Ny+y]);
		}
		fprintf(arq, "\n");
	}
}

int main(int argc, char **argv){

	int i, n, N, x, y, zm, I_inicial, I_0;
	int xp1, xm1, yp1, ym1;
	double t, a0, r, r_pa, Laplaciano, chance;
	double *S, *St;
	double *I, *It;
	double *R, *Rt;

	S= (double *) calloc(Nx*Ny, sizeof(double));
	I= (double *) calloc(Nx*Ny, sizeof(double));
	R= (double *) calloc(Nx*Ny, sizeof(double));
	St= (double *) calloc(Nx*Ny, sizeof(double));
	It= (double *) calloc(Nx*Ny, sizeof(double));
	Rt= (double *) calloc(Nx*Ny, sizeof(double));

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);

	/*for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			N= x*Ny + y;
			if(x == Nx/2 && y == Ny/2){
				S[N]= 0.0;
				I[N]= 1.0;
				R[N]= 0.0;
			}
			else{
				S[N]= 1.0;
				I[N]= 0.0;
				R[N]= 0.0;
			}
		}
	}*/
	I_0= 0;
	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			N= x*Ny + y;
			/*if(x == Nx/2 || x == Nx/2 + 1 || x == Nx/2 - 1 || x == Nx/2 + 2 || x == Nx/2 - 2){
				if(y == Ny/2 || y == Ny/2 + 1 || y == Ny/2 - 1 || y == Ny/2 + 2 || y == Ny/2 - 2){
					S[N]= 0.0;
					I[N]= 1.0;
					V[N]= 0.0;
				}
			}*/
			I_inicial= 10;
			chance= gsl_rng_uniform(w);
			if(chance < 0.0007){
				if(I_0 <= I_inicial){
					S[N]= 0.0;
					I[N]= 1.0;
					R[N]= 0.0;
					I_0++;
				}
			}
			else{
				S[N]= 1.0;
				I[N]= 0.0;
				R[N]= 0.0;
			}
		}
	}

	op(0, S, I, R);
	n= 1;
	r_pa= tf/NF;
	a0= r_pa;

	while(t < tf){
		for(x= 0; x< Nx; x++){ // LAPLACIANO ESPAÇO (primeiro passo RK)
			for(y= 0; y< Ny; y++){
				N= x*Ny + y;
				xp1= ((x+1)%Nx)*Ny+y;
				xm1= ((x-1+Nx)%Nx)*Ny+y;
				yp1= x*Ny+(y+1)%Ny;
				ym1= x*Ny+(y-1+Ny)%Ny;
				Laplaciano= S[xp1] + S[xm1] - 2.0*S[N] + S[yp1] + S[ym1] - 2.0*S[N];
				St[N]= S[N] + 0.5*dt*(D*Laplaciano - beta*S[N]*I[N]);
				Laplaciano= I[xp1] + I[xm1] - 2.0*I[N] +I[yp1] + I[ym1] - 2.0*I[N];
				It[N]= I[N] + 0.5*dt*(D*Laplaciano + beta*S[N]*I[N] - gamma*I[N]);
				Laplaciano= R[xp1] + R[xm1] - 2.0*R[N] + R[yp1] + R[ym1] - 2.0*R[N];
				Rt[N]= R[N] + 0.5*dt*(D*Laplaciano + gamma*I[N]);
			}
		}

		for(x= 0; x< Nx; x++){ // LAPLACIANO t (segundo passo RK)
			for(y= 0; y< Ny; y++){
				N= x*Ny + y;
				xp1= ((x+1)%Nx)*Ny+y;
				xm1= ((x-1+Nx)%Nx)*Ny+y;
				yp1= x*Ny+(y+1)%Ny;
				ym1= x*Ny+(y-1+Ny)%Ny;
				Laplaciano= St[xp1] + St[xm1] - 2.0*St[N] + St[yp1] + St[ym1] - 2.0*St[N];
				S[N]+= dt*(D*Laplaciano - beta*St[N]*It[N]);
				Laplaciano= It[xp1] + It[xm1] - 2.0*It[N] + It[yp1] + It[ym1] - 2.0*It[N];
				I[N]+= dt*(D*Laplaciano + beta*St[N]*It[N] - gamma*It[N]);
				Laplaciano= Rt[xp1] + Rt[xm1] - 2.0*Rt[N] + Rt[yp1] + Rt[ym1] - 2.0*Rt[N];
				R[N]+= dt*(D*Laplaciano + gamma*It[N]);
			}
		}

		t+= dt;
		if(t >= a0){
			op(n++, S, I, R);
			a0+= r_pa;
		}
	}
	free(S);
	free(I);
	free(R);
	free(St);
	free(It);
	free(Rt);
	gsl_rng_free(w);


	return 0;
}