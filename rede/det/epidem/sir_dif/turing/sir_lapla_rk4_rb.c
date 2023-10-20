// COMPILAÇÃO: gcc -O3 sir_lapla_rk4_rb.c -Wall -lm -lgsl -lgslcblas -o sir_lapla_rk4_rb
// EXECUÇÃO: ./sir_lapla_rk4_rb
// 1 = Suscetível, 2 = Infectado, 3 = Vacinado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define Nx 200  // tamanho da rede
#define Ny 200  // tamanho da rede
#define dt 0.1  // intervalo de tempo
#define D  1.0  // Parâmetro difusivo (mobilidade)
#define tf 5000 // tempo final
#define NF 5000 // número de arquivos de saída
#define beta 10.5 // 3.5    10.5
#define epsilon 0.3 // 0.1     0.3
#define d 0.3 // 0.1      0.3
#define mu 0.25 // 0.08       0.25
#define A 0.30 // 0.1        0.3
#define d11 1.8 // 0.6         1.8
#define d22 0.3 // 0.1			0.3
#define d33 1.8 // 0.1				0.3
#define d12 0.0 // 0.15				-0.0

void op(int n, double *S, double *I, double *V){
	int x, y;
	FILE *arq;
	char nome[100];
	sprintf(nome, "dat/sir-%d.dat", n);
	arq= fopen(nome, "w");

	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			fprintf(arq, "%e ", 1.0*S[x*Ny+y] + 2.0*I[x*Ny+y] + 3.0*V[x*Ny+y]);
		}
		fprintf(arq, "\n");
	}
	fclose(arq);
}

int main(int argc, char **argv){

	int I_inicial;
	int i, n, N, x, y, z;
	int xp1, xm1, yp1, ym1;
	double t, a0, r, r_pa, Laplaciano, Laplaciano_I;
	double chance;
	double *S, *St;
	double *I, *It;
	double *V, *Vt;

	S= (double *) calloc(Nx*Ny, sizeof(double));
	I= (double *) calloc(Nx*Ny, sizeof(double));
	V= (double *) calloc(Nx*Ny, sizeof(double));
	St= (double *) calloc(Nx*Ny, sizeof(double));
	It= (double *) calloc(Nx*Ny, sizeof(double));
	Vt= (double *) calloc(Nx*Ny, sizeof(double));

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);

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
			I_inicial= 20;
			chance= gsl_rng_uniform(w);
			if(chance < 0.0005){
				S[N]= 0.0;
				I[N]= 1.0;
				V[N]= 0.0;
			}
			else{
				S[N]= 1.0;
				I[N]= 0.0;
				V[N]= 0.0;
			}
		}
	}

	op(0, S, I, V);

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
				Laplaciano_I= I[xp1] + I[xm1] - 2.0*I[N] +I[yp1] + I[ym1] - 2.0*I[N];
				St[N]= S[N] + 0.5*dt*(d11*Laplaciano - beta*S[N]*I[N]*I[N] - d*S[N] + A + d12*Laplaciano_I);
				Laplaciano= I[xp1] + I[xm1] - 2.0*I[N] +I[yp1] + I[ym1] - 2.0*I[N];
				It[N]= I[N] + 0.5*dt*(d22*Laplaciano + beta*S[N]*I[N]*I[N] - epsilon*I[N] - d*I[N] - mu*I[N]);
				Laplaciano= V[xp1] + V[xm1] - 2.0*V[N] + V[yp1] + V[ym1] - 2.0*V[N];
				Vt[N]= V[N] + 0.5*dt*(d33*Laplaciano + epsilon*I[N] - d*V[N]);
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
				S[N]+= dt*(d11*Laplaciano - beta*S[N]*I[N]*I[N] - d*S[N] + A + d12*Laplaciano_I);
				Laplaciano= It[xp1] + It[xm1] - 2.0*It[N] + It[yp1] + It[ym1] - 2.0*It[N];
				I[N]+= dt*(d22*Laplaciano + beta*S[N]*I[N]*I[N] - epsilon*I[N] - d*I[N] - mu*I[N]);
				Laplaciano= Vt[xp1] + Vt[xm1] - 2.0*Vt[N] + Vt[yp1] + Vt[ym1] - 2.0*Vt[N];
				V[N]+= dt*(d33*Laplaciano + epsilon*I[N] - d*V[N]);
			}
		}

		t+= dt;
		/*if(t >= a0){
			op(n++, S, I, V);
			a0+= r_pa;
		}*/

		if(t >= a0){
			//if(n%50 == 0){
				op(n++, S, I, V);
				printf("%d\n", n);
			//}
			//else{
			//	n++;
			//}
			a0+= r_pa;
		}
	}
	free(S);
	free(I);
	free(V);
	free(St);
	free(It);
	free(Vt);
	gsl_rng_free(w);


	return 0;
}