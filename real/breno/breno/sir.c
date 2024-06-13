#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define NG 500	// gerações
#define NS 3 // espécies
#define I0 100
#define gamma 0.175
#define beta 1-gamma
/*#define PROBABILITY_OF_PREDATION \ // matriz de predação
const double pp[NS*NS]= {\
	0.0, 1.0, 0.0, \
	0.0, 0.0, 1.0, \
	1.0, 0.0, 0.0, \
};*/

int main(int argc, char **argv){
	if(argc != 6){
		printf("%s n_b.dat x_y.dat N M seed\n", argv[0]); // necessario o arquivo das vizinhanças, a posição dos indivíduos, número de indivíduos, vizinhanças médias e seed
		exit(1);
	}

	int i, j, k, t, ac, gd, nb, tp, n_b_max, ic, trig=0;
	int *p, *d, *n_b, *prey, S, I, R;
	//PROBABILITY_OF_PREDATION
	double at, act;
	double *x, *y;
	char name[100];
	FILE *file;

	int N = atoi(argv[3]);

	if(!(file = fopen(argv[1], "r"))){ // verifica se o arquivo n_b existe
		printf("cannot open file %s\n", argv[1]);
		exit(1);
	}
	n_b_max = 0;
	for(i = 0; i < N; i++){ // percorre os indivíduos
		if((fscanf(file, "%d", &nb) != 1)){ // verifica se o arquivo de dados é legível
			printf("1 - cannot read file %s\n", argv[1]);
			exit(1);
		}
		if(n_b_max < nb){
			n_b_max = nb; // pega o número de vizinhos do indivíduo que possui o maior número de vizinhos
		}
		for(j = 0; j < nb; j++){ // percorre os vizinhos dos indivíduos
			if((fscanf(file, "%d", &t) != 1)){
				printf("2 - cannot read file %s\n", argv[1]);
				exit(1);
			}
		}
	}
	fclose(file);
	n_b_max++;
	if((n_b = (int *) calloc(N*n_b_max, sizeof(int))) == NULL){ // aloca na memoria o número máximo de vizinhos permitidos nesse sistema
		printf("cannot allocate n_b\n");
		exit(1);
	}
	if(!(file = fopen(argv[1], "r"))){ // verifica de novo se o arquivo existe
		printf("cannot open file %s\n", argv[1]);
		exit(1);
	}
	for(i = 0; i < N; i++){ // 
		if((fscanf(file, "%d", &n_b[i*n_b_max]) != 1)){ // verifica se é possível ler o arquivo de dados. Por que i*N_b_max???
			printf("3 - cannot read file %s\n", argv[1]);
			exit(1);
		}
		k = n_b[i*n_b_max]; // não entendi o porquê de pegar i*n_b_max de n_b
		for(j = 1; j < k+1; j++){
			if((fscanf(file, "%d", &n_b[i*n_b_max+j]) != 1)){
				printf("4 - cannot read file %s\n", argv[1]);
				exit(1);
			}
		}
	}
	fclose(file);

	if((prey = (int *) calloc(n_b_max, sizeof(int))) == NULL){ // verifica se da pra alocar o vetor prey na memória
		printf("cannot allocate prey\n");
		exit(1);
	}

	gsl_rng_default_seed = atoi(argv[5]);
	gsl_rng *w = gsl_rng_alloc(gsl_rng_taus);

	if((p = (int *) calloc(N, sizeof(int))) == NULL){
		printf("cannot allocate p\n");
		exit(1);
	}

	if((d = (int *) calloc(NS, sizeof(int))) == NULL){
		printf("cannot allocate d\n");
		exit(1);
	}

	S= I= R= 0;
	for(i= 0; i< N; i++){
		p[i]= 1;
		S++;
	}
	for(i = 0; i < I0; i++){ // condição inicial. Para o modelo SIR, selecionar uns 10 indivíduos para começarem como infectados.
		ic= gsl_rng_uniform(w)*N;
		p[ic]= 2;
		S--;
		I++;
	}

	if(!(file = fopen(argv[2], "r"))){
		printf("cannot open file %s\n", argv[2]);
		exit(1);
	}
	if((x = (double *) calloc(N, sizeof(double))) == NULL){
		printf("cannot allocate x\n");
		exit(1);
	}
	if((y = (double *) calloc(N, sizeof(double))) == NULL){
		printf("cannot allocate y\n");
		exit(1);
	}
	for(i = 0; i < N; i++){
		if((fscanf(file, "%lf %lf", &x[i], &y[i]) != 2)){
			printf("cannot read file %s\n", argv[2]);
			exit(1);
		}
	}
	fclose(file);

	sprintf(name, "dat-%d/d-%d.dat", atoi(argv[4]), atoi(argv[5]));
	file = fopen(name, "w");
	fclose(file);

	for(t = 1; t <= NG; t++){
		gd = 0; // contador de passos de Monte Carlo
		while(gd < N){
			if(trig == 0){
				at = gsl_rng_uniform(w); // numero aleatorio
				ac = gsl_rng_uniform(w)*N; // indivíduo aleatório
				if(at < beta){ // prob de acontecer a mobilidade
					nb = gsl_rng_uniform(w)*(n_b[ac*n_b_max])+1; // pega um dos vizinhos do indivíduo ativo (individuo_ativo * número máximo de vizinhos) + 1
					nb = n_b[ac*n_b_max+nb]; // ?????? pega o vizinho
					if(p[ac] == 2){
						if(p[nb] == 1){
							p[nb]= p[ac];
							S--;
							I++;
						}
					}
					gd++;
				}else if(at >= beta){ // recuperação
					if(p[ac] == 2){
						p[ac]= 3;
						I--;
						R++;
					}
					gd++;
				}
			}
			else{
				gd= N;
			}
		}
		sprintf(name, "dat-%d/d-%d.dat", atoi(argv[4]), atoi(argv[5]));
		file = fopen(name, "a");
		fprintf(file, "%d %e %e %e\n", t, 1.0*S/N, 1.0*I/N, 1.0*R/N);
		fclose(file);
		if(I == 0){
			//t= NG+1;
			trig= 1;
		}
	}

	free(prey);
	free(n_b);
	free(x);
	free(y);
	free(p);
	gsl_rng_free(w);
	return 0;
}
