// Modelo de espalhamento de rumores estocástico SUVR POLARIZADO (Sensíveis, Indignados, Violentos e Relaxados) com presença de bolhas
// compilação: gcc suvr.c -O3 -Wall -lgsl -lgslcblas -lm -o suvr
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h>

#define N 10000  // número total de indivíduos (10000)
#define L 1.0   // tamanho do espaço
#define l 0.01 // tamanho do passo (0.01)
#define R 0.1 	 // tamanho da zona de interação - delimita o espaço de vizinhos
#define l_rep 0.01 // tamanho da distância de reprodução
#define NG 500 // número de gerações ###2 reduzindo o número para 400, só pra ficar mais rápido (1000)
#define pd 0.0  // probabilidade de morte (0.3)
#define pr 1.0  // probabilidade de reprodução (0.7)
#define M 140   // Limite de indivíduos do cluster
#define eta 1.0 // parametro direcional da mobilidade
#define l_inf 0.02 // raio de interação de infecção
#define l_polarizacao 0.03
#define I0 40 // número inicial de infectados #100
#define p_rec_u1 0.0 // probabilidade de se recuperar a partir do U        0.1
#define p_rec_v1 0.0 // probabilidade de se recuperar a partir do V	   0.05
#define p_rec_u2 0.0 // probabilidade de se recuperar a partir do U        0.1
#define p_rec_v2 0.0 // probabilidade de se recuperar a partir do V	   0.05
#define p_rec_v1u1 0.15 // probabilidade de re recuperar para U1, a partir de V1
#define p_rec_v2u2 0.15 // probabilidade de re recuperar para U2, a partir de V2
#define p_rec_rs 0.1    // probabilidade de re recuperar para S, a partir de R
#define p_u1su1 0.5 // probabilidade de gerar um U no contato entre S e U   0.5
#define p_u1u1v1 0.3 // probabilidade de gerar um V no contato entre U e U  0.3
#define p_u2su2 0.5 // probabilidade de gerar um U no contato entre S e U  0.5
#define p_u2u2v2 0.3 // probabilidade de gerar um V no contato entre U e U  0.3
#define p_v1su2 0.3  // probabilidade de gerar um U no contato entre V e S  0.3
#define p_v2su1 0.3  // probabilidade de gerar um U no contato entre V e S 0.3

// LEMBRETE: AJUSTAR AS PROBABILIDADES LÁ PRA BAIXO, NA PARTE DE INFECÇÃO

struct intble {
	int    s; // gênero
	int    c; // classe (morto= 0, sensíveis= 1, indignados1= 2, violentos1= 3, indignados2= 4, violentos2= 5, recuperadso= 6)
	double x; // posição em x
	double y; // posição em y
	double d; // distância euclidiana entre o indivíduo i e o j
};

void op(int t, struct intble *p){
	int i, j;
	FILE *arq;
	char clus[100];

	sprintf(clus, "dat/p-%d-%d.dat", 1, t); // Retorna uma string de acordo com a string de formatação Format
	arq= fopen(clus, "w"); // Abre o arquivo contido no primeiro parâmetro, usando o modo dado no segundo parâmetro. com o "w", o arquivo é criado se não existir.	
	for(i= 0; i< N; i++){
		if(p[i].s != 0){ // p[i].s acessa os elementos da variável s do struct
			fprintf(arq, "%e %e %d\n", p[i].x, p[i].y, p[i].c); 
		}
	}
	fclose(arq); 
}

int main(int argc, char **argv){
	int i, j, n, o, t, k, O, K, inf, i0, vio, i2; // k é o contador de indivíduos 
	int n_m= 0; // número de machos
	int n_f= 0; // número de femeas
	int n_s= 0; // número de suscetíveis
	int n_i= 0; // número de infectados
	double d, dx, dy, th, tht, x_cm, y_cm, x1, x2, y1, y2, t1, t2, D, DX, DY; // th é o ângulo theta
	double act, test, z, teste_inf, teste_i0, teste_rec, teste_inf1, teste_inf2, teste_i2; // test está relacionado com a mobilidade direcional, para obter um intervalo [-th, th], junto com z
	struct intble *p;
	FILE *arq;
	char clus[100];

	p= (struct intble *) calloc(N*N, sizeof(struct intble)); // Calloc aloca uma matriz na memória (número, tamanho)

	gsl_rng_default_seed= (argc == 2) ? atoi(argv[1]) : time(NULL);
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);// inicialização do gerador de número aleatório (algorítmo Taus)

	i0= 0;
	i2= 0;
	for(i= 0; i< N; i++){					// Atribui um gênero e posição para os individuos na rede
		p[i].s= 2.0*gsl_rng_uniform(w)+1.0; 
		p[i].x= gsl_rng_uniform(w);
		p[i].y= gsl_rng_uniform(w);
		p[i].d= L;
		if(I0 != i0){
			teste_i0= 1.3*gsl_rng_uniform(w);
			if(teste_i0 >= 1.0){
				p[i].c= 2;
				i0++;
			}
		}
		else if(I0 != i2){
			teste_i2= 1.3*gsl_rng_uniform(w);
			if(teste_i2 >= 1.0){
				p[i].c= 4;
				i2++;
			}
		}
		else{
			p[i].c= 1;
		}
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
			/*
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
			*/
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

				if(p[n].c == 2 || p[n].c == 3){
					for(j= 0; j< N; j++){
						if(p[j].c == 4 || p[j].c == 5){
							if(p[j].d < l_polarizacao){
								dx= p[j].x - p[n].x;
								dy= p[j].y - p[n].y;
								test = gsl_rng_uniform(w); 
								if(test < 0.5){
									z= 1;
								}
								else{
									z= -1;
								}
								tht= eta*z*M_PI*gsl_rng_uniform(w);
								p[n].x+= -dx + 0.4*l*cos(tht);
								if(p[n].x > L  ){ p[n].x-= L; }	 
								if(p[n].x < 0.0){ p[n].x+= L; }
								p[n].y+= -dy + 0.4*l*cos(tht);
								if(p[n].y > L  ){ p[n].y-= L; }
								if(p[n].y < 0.0){ p[n].y+= L; }
								j= N;
							}
						}
					}
				}
				if(p[n].c == 2 || p[n].c == 3){
					for(j= 0; j< N; j++){
						if(p[j].c == 2 || p[j].c == 3){
							if(p[j].d < l_polarizacao){
								dx= p[j].x - p[n].x;
								dy= p[j].y - p[n].y;
								test = gsl_rng_uniform(w); 
								if(test < 0.5){
									z= 1;
								}
								else{
									z= -1;
								}
								tht= eta*z*M_PI*gsl_rng_uniform(w);
								p[n].x+= +dx + 0.4*l*cos(tht);
								if(p[n].x > L  ){ p[n].x-= L; }	 
								if(p[n].x < 0.0){ p[n].x+= L; }
								p[n].y+= +dy + 0.4*l*cos(tht);
								if(p[n].y > L  ){ p[n].y-= L; }
								if(p[n].y < 0.0){ p[n].y+= L; }
								j= N;
							}
						}
					}
				}
				if(p[n].c == 4 || p[n].c == 5){
					for(j= 0; j< N; j++){
						if(p[j].c == 2 || p[j].c == 3){
							if(p[j].d < l_polarizacao){
								dx= p[j].x - p[n].x;
								dy= p[j].y - p[n].y;
								test = gsl_rng_uniform(w); 
								if(test < 0.5){
									z= 1;
								}
								else{
									z= -1;
								}
								tht= eta*z*M_PI*gsl_rng_uniform(w);
								p[n].x+= -dx + 0.4*l*cos(tht);
								if(p[n].x > L  ){ p[n].x-= L; }	 
								if(p[n].x < 0.0){ p[n].x+= L; }
								p[n].y+= -dy + 0.4*l*cos(tht);
								if(p[n].y > L  ){ p[n].y-= L; }
								if(p[n].y < 0.0){ p[n].y+= L; }
								j= N;
							}
						}
					}
				}
				if(p[n].c == 4 || p[n].c == 5){
					for(j= 0; j< N; j++){
						if(p[j].c == 4 || p[j].c == 5){
							if(p[j].d < l_polarizacao){
								dx= p[j].x - p[n].x;
								dy= p[j].y - p[n].y;
								test = gsl_rng_uniform(w); 
								if(test < 0.5){
									z= 1;
								}
								else{
									z= -1;
								}
								tht= eta*z*M_PI*gsl_rng_uniform(w);
								p[n].x+= +dx + 0.4*l*cos(tht);
								if(p[n].x > L  ){ p[n].x-= L; }	 
								if(p[n].x < 0.0){ p[n].x+= L; }
								p[n].y+= +dy + 0.4*l*cos(tht);
								if(p[n].y > L  ){ p[n].y-= L; }
								if(p[n].y < 0.0){ p[n].y+= L; }
								j= N;
							}
						}
					}
				}
				if(p[n].c == 1 || p[n].c == 6){
					test = gsl_rng_uniform(w); 
					if(test < 0.5){
						z= 1;
					}
					else{
						z= -1;
					}
					tht= eta*z*M_PI*gsl_rng_uniform(w); 
					p[n].x+= l*cos(tht);				 
					if(p[n].x > L  ){ p[n].x-= L; }	
					if(p[n].x < 0.0){ p[n].x+= L; }
					p[n].y+= l*sin(tht);				
					if(p[n].y > L  ){ p[n].y-= L; }
					if(p[n].y < 0.0){ p[n].y+= L; }
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

				inf= 0;
				vio= 0;
				if(p[n].c == 2 && p[n].s != 0){ // transm ativa
					for(j= 0; j< N; j++){
						teste_inf= gsl_rng_uniform(w);
						if(l_inf > p[j].d && teste_inf < p_u1su1 && p[j].c == 1 && p[j].s != 0 && inf != 5){ // Entre U1 e S gerando U1 || REVER ESSE p_usu
							inf++;				 
							p[j].c= 2;			  
						}
						if(l_inf > p[j].d && teste_inf < p_u1u1v1 && p[j].c == 2 && p[j].s != 0 && inf != 5){ // Entre U1 e U1 gerando V1
							p[j].c= 3;
							inf++;
						}
					}
				}
				if(p[n].c == 4 && p[n].s != 0){ 
					for(j= 0; j< N; j++){
						teste_inf= gsl_rng_uniform(w);
						if(l_inf > p[j].d && teste_inf < p_u2su2 && p[j].c == 1 && p[j].s != 0 && inf != 5){ // Entre U2 e S gerando U2 || REVER ESSE p_usu
							inf++;				 
							p[j].c= 4;			  
						}
						if(l_inf > p[j].d && teste_inf < p_u2u2v2 && p[j].c == 4 && p[j].s != 0 && inf != 5){ // Entre U2 e U2 gerando V2
							p[j].c= 5;
							inf++;
						}
					}
				}
				if(p[n].c == 3 && p[n].s != 0){
					for(j= 0; j< N; j++){
						teste_inf= gsl_rng_uniform(w);
						if(l_inf > p[j].d && teste_inf < p_v1su2 && p[j].c == 1 && p[j].s != 0 && vio != 4){
							vio++;			 
							p[j].c= 4;
						}
					}
				}
				if(p[n].c == 5 && p[n].s != 0){
					for(j= 0; j< N; j++){
						teste_inf= gsl_rng_uniform(w);
						if(l_inf > p[j].d && teste_inf < p_v2su1 && p[j].c == 1 && p[j].s != 0 && vio != 4){
							vio++;			 
							p[j].c= 2;
						}
					}
				}
				printf("%d com geração igual a %d\n", K, t);
				printf("Pŕoximo indivíduo. \n");
				// teste fim

				for(j= 0; j< N; j++){ // esta é a distância de reprodução ou l_rep
					if(l_rep >= p[o].d && K <= M){ 
						if(p[j].s == 0){ 	// ****** testei com (&& k <= M && l >= d)		// Se ouver algum espaço livre na rede (indiv. < N), executa a rep.
							p[j].s= 2.0*gsl_rng_uniform(w)+1.0; // Esse bloquinho escolhe o gênero dos novos individ.
							p[j].c= 1;
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
							}else{      	  					// Partindo de uma fêmea também
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
// recuperação - começo
				if(p[n].c == 2 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_u1){
						p[n].c= 6;
					}
				}
				if(p[n].c == 3 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_v1){
						p[n].c= 6;
					}
				}
				if(p[n].c == 4 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_u2){
						p[n].c= 6;
					}
				}
				if(p[n].c == 5 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_v2){
						p[n].c= 6;
					}
				}
				if(p[n].c == 3 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_v1u1){
						p[n].c= 2;
					}
				}
				if(p[n].c == 5 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_v2u2){
						p[n].c= 4;
					}
				}
				if(p[n].c == 6 && p[n].s != 0){
					teste_rec= gsl_rng_uniform(w);
					if(teste_rec < p_rec_rs){
						p[n].c= 1;
					}
				}
// recuperação - fim
			}else{
// morte - começo
				//E se eu colocar uma probabilidade disso acontecer, dentro da probabilidade entre rep e morte?
				//act;
				//if (act < pd2){
				if(p[n].s == 1){n_f--;} // redução no número de indivíduos F ou M
				else{n_m--;}
				p[n].s= 0;
				p[n].c= 0; // liberação de um espaço
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
