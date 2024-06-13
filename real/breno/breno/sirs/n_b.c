// gera a vizinhanca dos individuos
// gcc nb.c -O3 -Wall -lgsl -lgslcblas -lm -o nb
// ./nb densidade_vizinhanda numero_de_individuos seed

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_rng.h>

#define L 1.0
#define l_i 0.01

struct intble{
	int nb;
	double x;
	double y;
	double d;
};

int main(int argc, char **argv){
    int i, j, N, cont, nb, o;
    double dx, dy, dd;
	struct intble *p;
    FILE *arq;
    char name[100];
	char name2[100];

	nb= atoi(argv[1]);
    N= atoi(argv[2]);
	p= (struct intble *) calloc(N*N*N, sizeof(struct intble));


    gsl_rng_default_seed = atoi(argv[3]);
	gsl_rng *w = gsl_rng_alloc(gsl_rng_taus);

    sprintf(name, "xy/xy-%d-%ld.dat", N, gsl_rng_default_seed);
	arq = fopen(name, "w");
    for(i= 0; i< N; i++){
        p[i].x= gsl_rng_uniform(w);
        p[i].y= gsl_rng_uniform(w);
		p[i].d= 1.0;
		p[i].nb= 0;
        fprintf(arq, "%e %e\n", p[i].x, p[i].y);
    }
    fclose(arq);

	sprintf(name2, "nb/nb-%d-%d-%ld.dat", nb, N, gsl_rng_default_seed);
	arq = fopen(name2, "w");

	o= 0; // 
	dd= L;
	for(i= 0; i< N; i++){
		nb= 0;
		for(j= 0; j< N; j++){
			p[j].d = L;
			if(p[j].d != 0){
				dx= fabs(p[j].x - p[i].x);		
				if(dx > 0.5*L){ dx-= L; } 		
				dy= fabs(p[j].y - p[i].y); 		
				if(dy > 0.5*L){ dy-= L; } 		
				if(sqrt(dx*dx+dy*dy) < dd){
					dd= sqrt(dx*dx+dy*dy); 		
					o= j;
					p[j].d= dd;
				}
			}
			p[j].nb= 0;
			if(p[j].d< l_i && p[j].d != 0){
				nb++;
				p[j].nb= j;
			}
		}
		fprintf(arq, "%d ", nb);
		for(j= 0; j< N; j++){
			if(p[j].nb != 0){
				fprintf(arq, "%d ", p[j].nb);
			}
		}
		fprintf(arq, "\n");
	}
	fclose(arq);
	
	free(p);
    gsl_rng_free(w);
    return 0;
}