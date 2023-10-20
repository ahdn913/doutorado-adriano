// COMPILAÇÃO: gcc -O3 rps.c -Wall -lm -lgsl -lgslcblas -o rps
// EXECUÇÃO: ./rps
// Espécies A(1) -> B(2) -> C(3)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define Nx 100  // tamanho da rede
#define Ny 100  // tamanho da rede
#define dt 0.1  // intervalo de tempo
#define D  0.3  // Parâmetro difusivo (mobilidade)
#define tf 300 // tempo final
#define NF 3000 // número de arquivos de saída
#define beta 0.85
#define gamma 0.15
#define Ne 3 // número de espécies

void op(int n, double *A, double *B, double *C){ // printa a matriz de estado de cada geração
	int x, y;
	FILE *arq;
	char nome[100];
	sprintf(nome, "dat_plot/sir-%d.dat", n);
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

	int i, n, N, x, y, zm, inicial, cont_A, cont_B, cont_C;
	int xp1, xm1, yp1, ym1;
	double t, a0, r, r_pa, Laplaciano, chance;
	double *A, *At;
	double *B, *Bt;
	double *C, *Ct;
	FILE *arquivo;
	char nome1[100];

	A= (double *) calloc(Nx*Ny, sizeof(double));
	B= (double *) calloc(Nx*Ny, sizeof(double));
	C= (double *) calloc(Nx*Ny, sizeof(double));
	At= (double *) calloc(Nx*Ny, sizeof(double));
	Bt= (double *) calloc(Nx*Ny, sizeof(double));
	Ct= (double *) calloc(Nx*Ny, sizeof(double));

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);

	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			N= x*Ny + y;
			inicial= gsl_rng_uniform(w)*3;
			switch(Ne){
				case 0:
					A[N]= 1.0;
					B[N]= 0.0;
					C[N]= 0.0;
				break;
				case 1:
					A[N]= 0.0;
					B[N]= 1.0;
					C[N]= 0.0;
				break;
				default:
					A[N]= 0.0;
					B[N]= 0.0;
					C[N]= 1.0;
				break;
			}
		}
	}

	t= 0.0;
	n= 0;
	sprintf(nome1, "dat_curve/sir-%d.dat", n);
	arquivo= fopen(nome1, "w");
	cont_A= cont_B= cont_C= 0;
	// for(x= 0; x< Nx; x++){
	// 	for(y= 0; y< Ny; y++){
	// 		if(A[N] >= 1.0 && A[N] < 2.0){
	// 			cont_A++;
	// 		}
	// 		else if(B[N] >= 2.0 && B[N] < 3.0){
	// 			cont_B++;
	// 		}
	// 		else if(C[N] >= 3.0 && C[N] < 4.0){
	// 			cont_C++;
	// 		}
	// 	}
	// }
	for(x= 0; x< Nx; x++){
		for(y= 0; y< Ny; y++){
			if(A[N] != 0.0){
				cont_A++;
			}
			else if(B[N] != 0.0){
				cont_B++;
			}
			else if(C[N] != 0.0){
				cont_C++;
			}
		}
	}
	fprintf(arquivo, "%e %d %d %d\n", t, cont_A, cont_B, cont_C);
	fclose(arquivo);
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
				At[N]= A[N] + 0.5*dt*(D*Laplaciano - A[N]*A[N] + A[N] - beta*C[N]*A[N]); // basta dividir por (Nx*Ny)
				Laplaciano= B[xp1] + B[xm1] - 2.0*B[N] +B[yp1] + B[ym1] - 2.0*B[N];
				Bt[N]= B[N] + 0.5*dt*(D*Laplaciano -B[N]*B[N] + B[N] - beta*B[N]*A[N]);
				Laplaciano= C[xp1] + C[xm1] - 2.0*C[N] + C[yp1] + C[ym1] - 2.0*C[N];
				Ct[N]= C[N] + 0.5*dt*(D*Laplaciano - C[N]*C[N] + C[N] - beta*C[N]*B[N]);
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
				A[N]+= dt*(D*Laplaciano - At[N]*At[N] + At[N] - beta*Ct[N]*At[N]);
				Laplaciano= Bt[xp1] + Bt[xm1] - 2.0*Bt[N] + Bt[yp1] + Bt[ym1] - 2.0*Bt[N];
				B[N]+= dt*(D*Laplaciano -Bt[N]*Bt[N] + Bt[N] - beta*Bt[N]*At[N]);
				Laplaciano= Ct[xp1] + Ct[xm1] - 2.0*Ct[N] + Ct[yp1] + Ct[ym1] - 2.0*Ct[N];
				C[N]+= dt*(D*Laplaciano - Ct[N]*Ct[N] + Ct[N] - beta*Ct[N]*Bt[N]);
			}
		}

		t+= dt;
		if(t >= a0){
			op(n++, A, B, C);
			a0+= r_pa;
		}

		t= 0.0;
		sprintf(nome1, "dat_curve/sir-%d.dat", n);
		arquivo= fopen(nome1, "a");
		cont_A= cont_B= cont_C= 0;
		// for(x= 0; x< Nx; x++){
		// 	for(y= 0; y< Ny; y++){
		// 		if(A[N] >= 1.0 && A[N] < 2.0){
		// 			cont_A++;
		// 		}
		// 		else if(B[N] >= 2.0 && B[N] < 3.0){
		// 			cont_B++;
		// 		}
		// 		else if(C[N] >= 3.0 && C[N] < 4.0){
		// 			cont_C++;
		// 		}
		// 	}
		// }
		for(x= 0; x< Nx; x++){
			for(y= 0; y< Ny; y++){
				if(A[N] != 0.0){
					cont_A++;
				}
				else if(B[N] != 0.0){
					cont_B++;
				}
				else if(C[N] != 0.0){
					cont_C++;
				}
			}
		}
		fprintf(arquivo, "%e %d %d %d\n", t, cont_A, cont_B, cont_C);
		fclose(arquivo);

	}
	free(A);
	free(B);
	free(C);
	free(At);
	free(Bt);
	free(Ct);
	gsl_rng_free(w);
	return 0;
}