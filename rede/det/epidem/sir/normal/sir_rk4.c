// COMPILAÇÃO: gcc sir_rk4.c -Wall -lm -o sir_rk4
// EXECUÇÃO: ./sir_rk4 s0 i0 r0 dt tf np (suscetíveis, infectados e vacinados iniciais, intervalo de tempo, tempo total e número de pontos)
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

double f(double S, double I){ // primeira eq dif S
	double beta= 0.85;
	return (-beta*S*I);
}

double g(double S, double I){ // segunda eq dif I
	double beta= 0.85;
	double gamma= 0.15;
	return (beta*S*I - gamma*I);
}

double h(double I){ // terceira eq dif R
	double gamma= 0.15;
	return (gamma*I);
}

int main(int argc, char **argv){
	if (argc != 7){
		printf("Erro, dá uma lida na execução do programa, número errado de entradas.\n");
		exit(1);
	};

	int i, j, k, n, np;
	double S, I, R, t, tf, dt, k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4;

	S= atof(argv[1]);
	I= atof(argv[2]);
	R= atof(argv[3]);
	dt= atof(argv[4]);
	tf= atof(argv[5]);
	np= atoi(argv[6]);
	n=tf/(np*dt);
	FILE *arq;
	arq= fopen("sir_rk4.dat", "w");
	t= 0.0;

	for(i= 0; i< np; i++){
		for(j= 0; j< n; j++){ 
			k1= f(S, I); // inclinação da reta tangente
			l1= g(S, I); // inclinação da reta tangente
			m1= h(I);
			k2= f(S+0.5*dt*k1, I+0.5*dt*l1);
			l2= g(S+0.5*dt*k1, I+0.5*dt*l1);
			m2= h(I+0.5*dt*k1);
			k3= f(S+0.5*dt*k2, I+0.5*dt*l2);
			l3= g(S+0.5*dt*k2, I+0.5*dt*l2);
			m3= h(I+0.5*dt*k2);
			k4= f(S+dt*k3, I+dt*l3);
			l4= g(S+dt*k3, I+dt*l3);
			m4= h(I+dt*k3);
			S+= dt*(k1+2.0*k2+2.0*k3+k4)/6.0; //passo temporal
			I+= dt*(l1+2.0*l2+2.0*l3+l4)/6.0; //passo temporal
			R+= dt*(m1+2.0*m2+2.0*m3+m4)/6.0;
			t+= dt;
		}
		fprintf(arq, "%e %e %e %e\n", t, S, I, R);
	}

	fclose(arq);
	return 0;
}
