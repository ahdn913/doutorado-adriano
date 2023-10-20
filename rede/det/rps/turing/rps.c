// COMPILAÇÃO: gcc -O3 rps.c -Wall -lm -lgsl -lgslcblas -o rps
// EXECUÇÃO: ./rps
// Espécies A(1) -> B(2) -> C(3)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define Nx 300  // tamanho da rede
#define Ny 300  // tamanho da rede
#define dt 0.1  // intervalo de tempo
#define d1 1.2 // Parâmetro difusivo (mobilidade)
#define d2 1.2  // Parâmetro difusivo (mobilidade)
#define d3 1.2  // Parâmetro difusivo (mobilidade)
#define d4 0.1  // Parâmetro difusivo (mobilidade)
#define tf 2000 // tempo final
#define NF 2000 // número de arquivos de saída
#define beta11 0.0 // predação
#define beta22 0.0 // predação
#define beta33 0.0 // predação
#define beta21 1.0 // predação
#define beta12 0.0 // predação
#define beta23 0.0 // predação
#define beta31 0.0 // predação
#define beta32 1.0 // predação
#define beta13 1.0 // predação
#define r1 2.0 // reprodução (associado ao espaço vazio)
#define r2 2.0 // reprodução (associado ao espaço vazio)
#define r3 2.0 // reprodução (associado ao espaço vazio)
#define Ne 3 	  // número de espécies (incluindo o espaço vazio)

void op(int n, double *A, double *B, double *C){ // printa a matriz de estado de cada geração
	int x, y;
	FILE *arq;
	char nome[100];
	sprintf(nome, "dat_plot/rps-%d.dat", n);
	arq= fopen(nome, "w");

	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			fprintf(arq, "%e ", 1.0*A[x*Ny+y] + 2.0*B[x*Ny+y] + 3.0*C[x*Ny+y]);
		}
		fprintf(arq, "\n");
	}
	fclose(arq);
}

int main(int argc, char **argv){

	int i, n, N, x, y, zm, inicial, cont_a, cont_b, cont_c;
	int xp1, xm1, yp1, ym1;
	double t, a0, r_pa, Laplaciano, chance;
	double *A, *At;
	double *B, *Bt;
	double *C, *Ct;
	// double *O, *Ot;
	FILE *arquivo;
	char nome1[100];

	// como fazer isso aqui embaixo ser generalizado????
	// Se eu associar números diferentes para uma mesma variável, dá certo?

	A= (double *) calloc(Nx*Ny, sizeof(double));
	B= (double *) calloc(Nx*Ny, sizeof(double));
	C= (double *) calloc(Nx*Ny, sizeof(double));
	// O= (double *) calloc(Nx*Ny, sizeof(double));
	At= (double *) calloc(Nx*Ny, sizeof(double));
	Bt= (double *) calloc(Nx*Ny, sizeof(double));
	Ct= (double *) calloc(Nx*Ny, sizeof(double));
	// Ot= (double *) calloc(Nx*Ny, sizeof(double));

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);

	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			N= x*Ny + y;
			inicial= gsl_rng_uniform(w)*Ne;
			switch(inicial){
				case 0:
					A[N]= 1.0;
					B[N]= 0.0;
					C[N]= 0.0;
					// O[N]= 0.0;
				break;
				case 1:
					A[N]= 0.0;
					B[N]= 1.0;
					C[N]= 0.0;
					// O[N]= 0.0;
				break;
				default:
					A[N]= 0.0;
					B[N]= 0.0;
					C[N]= 1.0;
					// O[N]= 1.0;
				break;
			}
		}
	}

	gsl_rng_free(w);
	t= 0.0;
	n= 0;
	op(0, A, B, C);
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
				Laplaciano= A[xp1] + A[xm1] - 2.0*A[N] + A[yp1] + A[ym1] - 2.0*A[N];
				At[N]= A[N] + 0.5*dt*(d1*Laplaciano + r1*A[N]*(1.0 - A[x*Ny+y] - B[x*Ny+y] - C[x*Ny+y]) - beta11*A[N]*A[N] - beta12*A[N]*B[N] - beta13*A[N]*C[N]); // basta dividir por (Nx*Ny)
				Laplaciano= B[xp1] + B[xm1] - 2.0*B[N] +B[yp1] + B[ym1] - 2.0*B[N];
				Bt[N]= B[N] + 0.5*dt*(d2*Laplaciano + r2*B[N]*(1.0 - A[x*Ny+y] - B[x*Ny+y] - C[x*Ny+y]) - beta21*B[N]*A[N] - beta22*B[N]*B[N] - beta23*B[N]*C[N]);
				Laplaciano= C[xp1] + C[xm1] - 2.0*C[N] + C[yp1] + C[ym1] - 2.0*C[N];
				Ct[N]= C[N] + 0.5*dt*(d3*Laplaciano + r3*C[N]*(1.0 - A[x*Ny+y] - B[x*Ny+y] - C[x*Ny+y]) - beta31*C[N]*A[N] - beta32*C[N]*B[N] - beta33*C[N]*C[N]);
				// Laplaciano= O[xp1] + O[xm1] - 2.0*O[N] + O[yp1] + O[ym1] - 2.0*O[N];
				// Ot[N]= O[N] + 0.5*dt*(D*Laplaciano - r*O[N]*(A[N] + B[N] + C[N]) + beta*(A[N]*B[N] + B[N]*C[N] + C[N]*A[N]));
			}
		}

		for(x= 0; x< Nx; x++){ // LAPLACIANO t (segundo passo RK)
			for(y= 0; y< Ny; y++){
				N= x*Ny + y;
				xp1= ((x+1)%Nx)*Ny+y;
				xm1= ((x-1+Nx)%Nx)*Ny+y;
				yp1= x*Ny+(y+1)%Ny;
				ym1= x*Ny+(y-1+Ny)%Ny;
				Laplaciano= At[xp1] + At[xm1] - 2.0*At[N] + At[yp1] + At[ym1] - 2.0*At[N];
				A[N]+= dt*(d1*Laplaciano + r1*At[N]*(1.0 - At[x*Ny+y] - Bt[x*Ny+y] - Ct[x*Ny+y]) - beta11*At[N]*At[N] - beta12*At[N]*Bt[N] - beta13*At[N]*Ct[N]);
				Laplaciano= Bt[xp1] + Bt[xm1] - 2.0*Bt[N] + Bt[yp1] + Bt[ym1] - 2.0*Bt[N];
				B[N]+= dt*(d2*Laplaciano + r2*Bt[N]*(1.0 - At[x*Ny+y] - Bt[x*Ny+y] - Ct[x*Ny+y]) - beta21*Bt[N]*At[N] - beta22*Bt[N]*Bt[N] - beta23*Bt[N]*Ct[N]);
				Laplaciano= Ct[xp1] + Ct[xm1] - 2.0*Ct[N] + Ct[yp1] + Ct[ym1] - 2.0*Ct[N];
				C[N]+= dt*(d3*Laplaciano + r3*Ct[N]*(1.0 - At[x*Ny+y] - Bt[x*Ny+y] - Ct[x*Ny+y]) - beta31*Ct[N]*At[N] - beta32*Ct[N]*Bt[N] - beta33*Ct[N]*Ct[N]);
				// Laplaciano= Ot[xp1] + Ot[xm1] - 2.0*Ot[N] + Ot[yp1] + Ot[ym1] - 2.0*Ot[N];
				// O[N]+= dt*(D*Laplaciano - r*Ot[N]*(At[N] + Bt[N] + Ct[N]) + beta*(At[N]*Bt[N] + Bt[N]*Ct[N] + Ct[N]*At[N]));
			}
		}

		// t+= dt;
		// if(t >= a0){
		// 	op(n++, A, B, C, O);
		// 	a0+= r_pa;
		// }

		t+= dt;
		if(t >= a0){
			op(n++, A, B, C);
			a0+= r_pa;
		}
	}
	free(A);
	free(B);
	free(C);
	// free(O);
	free(At);
	free(Bt);
	free(Ct);
	// free(Ot);
	return 0;
}