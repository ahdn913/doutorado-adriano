// COMPILAÇÃO: gcc siv_rk4.c -Wall -lm -o siv_rk4
// EXECUÇÃO: ./siv_rk4 s0 i0 v0 dt tf np (suscetíveis, infectados e vacinados iniciais, intervalo de tempo, tempo total e número de pontos)
// 1 = Suscetível, 2 = Infectado, 4 = Vacinado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

double f(double S, double I){ // primeira eq dif S
	double beta= 0.12; //0.2
	double epsilon= 0.008; //0.01
	return (-beta*S*I - epsilon*S);
}

double g(double S, double I){ // segunda eq dif I
	double beta= 0.12; // 0.2
	return (beta*S*I);
}

double h(double S){ // terceira eq dif V
	double epsilon= 0.008; // 0.01
	return (epsilon*S);
}

int main(int argc, char **argv){
	if (argc != 7){
		printf("Erro, da uma lida na execução do programa, número errado de entradas.\n");
		exit(1);
	};

	int i, j, k, n, np;
	double S, I, V, t, tf, dt, k1, k2, k3, k4, l1, l2, l3, l4, m1, m2, m3, m4;

	S= atof(argv[1]);
	I= atof(argv[2]);
	V= atof(argv[3]);
	dt= atof(argv[4]);
	tf= atof(argv[5]);
	np= atoi(argv[6]);
	n=tf/(np*dt);
	FILE *arq;
	arq= fopen("siv_rk4.dat", "w");
	t= 0.0;

	for(i= 0; i< np; i++){
		for(j= 0; j< n; j++){ 
			k1= f(S, I); // inclinação da reta tangente
			l1= g(S, I); // inclinação da reta tangente
			m1= h(S);
			k2= f(S+0.5*dt*k1, I+0.5*dt*l1);
			l2= g(S+0.5*dt*k1, I+0.5*dt*l1);
			m2= h(S+0.5*dt*k1);
			k3= f(S+0.5*dt*k2, I+0.5*dt*l2);
			l3= g(S+0.5*dt*k2, I+0.5*dt*l2);
			m3= h(S+0.5*dt*k2);
			k4= f(S+dt*k3, I+dt*l3);
			l4= g(S+dt*k3, I+dt*l3);
			m4= h(S+dt*k3);
			S+= dt*(k1+2.0*k2+2.0*k3+k4)/6.0; //passo temporal
			I+= dt*(l1+2.0*l2+2.0*l3+l4)/6.0; //passo temporal
			V+= dt*(m1+2.0*m2+2.0*m3+m4)/6.0;
			t+= dt;
		}
		fprintf(arq, "%e %e %e %e\n", t, S, I, V);
	}

	fclose(arq);
	return 0;
}
