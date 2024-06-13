#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define N 10000  // número total de indivíduos (10000)
#define L 1.0   // tamanho do espaço
#define l 0.01 // tamanho do passo (0.01)
#define R 0.1 	 // tamanho da zona de interação
#define l_rep 0.01 // tamanho da distância de reprodução
#define NG 500 // número de gerações #2 reduzindo o número para 400, só pra ficar mais rápido (1000)
#define pd 0.3  // probabilidade de morte (0.3)
#define pr 0.7  // probabilidade de reprodução (0.7)
#define M 160   // Limite de indivíduos do cluster
#define eta 1.0 // parametro direcional da mobilidade

struct intble {
	int    s; // gênero
	double x; // posição em x
	double y; // posição em y
	double d; // 
};

void op(int t, struct intble *p){
	int i, j;
	FILE *arq;
	char clus[100];

	sprintf(clus, "dat/p-%d-%d.dat", 1, t); // Retorna uma string de acordo com a string de formatação Format
	arq= fopen(clus, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
	for(i= 0; i< N; i++){
		if(p[i].s != 0){ // p[i].s acessa os elementos da variável s do struct
			fprintf(arq, "%e %e %d\n", p[i].x, p[i].y, p[i].s); 
		}
	}
	fclose(arq); 
}

int main(int argc, char **argv){
	int i, j, n, o, t, k, O, K; // k é o contador de indivíduos 
	int n_m= 0; // número de machos
	int n_f= 0; // número de femeas
	double d, dx, dy, th, tht, x_cm, y_cm, x1, x2, y1, y2, t1, t2, D, DX, DY; // th é o ângulo theta
	double act, test, z; // test está relacionado com a mobilidade direcional, para obter um intervalo [-th, th], junto com z
	struct intble *p;
	FILE *arq;
	char clus[100];

	p= (struct intble *) calloc(N*N, sizeof(struct intble)); // Calloc aloca uma matriz na memória (número, tamanho)

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);// inicialização do gerador de número aleatório (algorítmo Taus)

	for(i= 0; i< N; i++){					// Atribui um gênero e posição para os individuos na rede
		p[i].s= 2.0*gsl_rng_uniform(w)+1.0; 
		p[i].x= gsl_rng_uniform(w);
		p[i].y= gsl_rng_uniform(w);
		p[i].d= L;
	}

	op(0, p); 				// Escreve a condição inicial dos indivíduos

	for(i= 0; i< N; i++){ 	// adiciona os indivíduos femininos e masculinos
		if(p[i].s == 1){
			n_f++;
		}else{
			n_m++;
		}
	}
	arq= fopen("dst.dat", "a");
	fprintf(arq, "%e %e\n", 1.0*n_f/N, 1.0*n_m/N); 
	fclose(arq);
	
	for(t= 1; t< NG; t++){ 		// Do tempo inicial até o número de gerações
		for(i= 0; i< N; i++){ 	// Do primeiro indivíduo até o último da população
			do{
				n= gsl_rng_uniform(w)*N; //
			}while(p[n].s == 0);
// mobilidade - começo
			test = gsl_rng_uniform(w); // este bloco de if e else serve para determinar z, que influi o range de theta [-tht, tht]
			if(test < 0.5){
				z= 1;
			}
			else{
				z= -1;
			}
			tht= eta*z*M_PI*gsl_rng_uniform(w); // Ângulo em radianos *********************** MUDEI TH PARA THT
			p[n].x+= l*cos(tht);				 // Tamanho do passo dado no eixo x
			if(p[n].x > L  ){ p[n].x-= L; }	 // Aqui estão listadas as condições de borda periódica
			if(p[n].x < 0.0){ p[n].x+= L; }
			p[n].y+= l*sin(tht);				 // É realizada a mesma coisa pro eixo y
			if(p[n].y > L  ){ p[n].y-= L; }
			if(p[n].y < 0.0){ p[n].y+= L; }
// mobilidade - fim
// ação
			act= gsl_rng_uniform(w); 	// gera-se um número aleatório entre 0 e 1
			if(act < pr){ 				// se este número for menor que a prob de repr., então têm-se a reprodução
// reprodução - começo
				o= 0; // 
				d= L;
				k=0;
				for(j= 0; j< N; j++){
					p[j].d = L;
					if(p[j].s != p[n].s && p[j].s != 0){// se o sexo do indivíduo escolhido for diferente do outro ind.
						dx= fabs(p[j].x - p[n].x);		// fabs te dá o módulo de um float/double
						if(dx > 0.5*L){ dx-= L; } 		// condição de borda em x
						dy= fabs(p[j].y - p[n].y); 		// diferença na distância no eixo y entre um indiv. e outro.
						if(dy > 0.5*L){ dy-= L; } 		// condição de borda em y
						if(sqrt(dx*dx+dy*dy) < d){
							d= sqrt(dx*dx+dy*dy); 		// distância do segmento r entre um indívíduo e outro
							o= j;
							//if(l >= d){				// ****** Testei com if(l >= d){k++;} mas parece não respeitar o M
							k++;				// O problema está na contagem de k. Fiz um printf e o valor max foi 6
							p[j].d= d;
							//printf("%d \n", k); // Ao retirar essa condição if, o valor max vai para 15
							//} // 
						}
					}
				}

				// teste inicio
				O= 0;
				D= L;
				K= 0;
				for(j= 0; j< N; j++){
					DX= fabs(p[j].x - p[n].x); //p[n]
					if(DX > 0.5*L){ DX-= L; }
					DY= fabs(p[j].y - p[n].y); //p[n]
					if(DY > 0.5*L){ DY-= L; }
					if(sqrt(DX*DX+DY*DY) < L){
						D= sqrt(DX*DX+DY*DY);
						O= j;
						if(D <= R && p[j].s != 0){   		// zona de interação R 
							K++;
						}
						//printf("%d\n", K);
					}
				}
				printf("%d com geração igual a %d\n", K, t);
				printf("Pŕoximo indivíduo. \n");
				// teste fim

				for(j= 0; j< N; j++){ // esta é a distância de reprodução ou l_rep
					if(l_rep >= p[o].d && K <= M){ 
						if(p[j].s == 0){ 	// ****** testei com (&& k <= M && l >= d)		// Se ouver algum espaço livre na rede (indiv. < N), executa a rep.
							p[j].s= 2.0*gsl_rng_uniform(w)+1.0; // Esse bloquinho escolhe o gênero dos novos individ.
							if(p[j].s == 1){
								n_f++;
							}else{
								n_m++;
							}
							th= 2.0*M_PI*gsl_rng_uniform(w);
							if(p[o].s == 1){  					// p[o]     Aqui se dá a reprodução, partindo de uma fêmea
								p[j].x= p[o].x + l*cos(th);  	// p[o]     posição do novo indivíduo
								if(p[j].x > L  ){ p[j].x-= L; } // condição de borda periódica
								if(p[j].x < 0.0){ p[j].x+= L; }
								p[j].y= p[o].y + l*sin(th);		// p[o]
								if(p[j].y > L  ){ p[j].y-= L; }
								if(p[j].y < 0.0){ p[j].y+= L; }
							}else{      	  					// Partindo de um macho
								p[j].x= p[n].x + l*cos(th);		// p[n] Aqui, utilizamos p[n] porque os individuos devem nascer
								if(p[j].x > L  ){ p[j].x-= L; } // próximos da mãe. Então se p[o].s != 1, é um macho.
								if(p[j].x < 0.0){ p[j].x+= L; }
								p[j].y= p[n].y + l*sin(th);		// p[n]
								if(p[j].y > L  ){ p[j].y-= L; }
								if(p[j].y < 0.0){ p[j].y+= L; }
							}
							j= N; // Quando n interage com j, o loop acaba e passa a vez pro próximo indivíduo
						}
					}
				}
// reprodução - fim
			}else{
// morte - começo
				//E se eu colocar uma probabilidade disso acontecer, dentro da probabilidade entre rep e morte?
				//act;
				//if (act < pd2){
				if(p[n].s == 1){n_f--;} // redução no número de indivíduos F ou M
				else{n_m--;}
				p[n].s= 0; // liberação de um espaço
				//}
// morte - fim
			}
		}
		op(t, p);
		printf("%d/%d\n", t, NG);
		arq= fopen("dst.dat", "a");
		fprintf(arq, "%e %e\n", 1.0*n_f/N, 1.0*n_m/N);
		fclose(arq);
	}
	gsl_rng_free(w);
	free(p);
	return 0;
}
