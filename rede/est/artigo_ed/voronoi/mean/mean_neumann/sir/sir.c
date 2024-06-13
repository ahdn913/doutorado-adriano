// gcc sir.c -O3 -Wall -lgsl -lgslcblas -lm -o sir
// ./sir
// for I in $(seq 1 100); do echo $I; ./sir $I ;done
// 1 = Suscetível, 2 = Infectado, 3 = Recuperado, 4 = Vacinado

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gsl/gsl_rng.h> 

#define L 100 // L = Ni = Nj
#define Ni 100   // tamanho da rede em i
#define Nj 100   // tamanho da rede em j
#define NG 1000
#define gamma 0.175
#define beta 1-gamma

int phi[L*L];

void op(int t){
	int i, j; 
	double x, y;
	FILE *arq;
	char nome[100];

	sprintf(nome, "plot_dat/sir-%d.dat", t);
	arq= fopen(nome, "w");
	for(j= 0; j< Nj; j++){
		for(i= 0; j< Ni; i++){

/*// Rede honeycomb 
 		int desl;
 		if(j%2==0){
			desl= 1;
		}else{
			desl= 0;
		}		

		for(i= 0; i< Ni; i++){

// Rede honeycomb			

			x= 1.0*i+desl;
			y= 1.0*j*(sqrt(3.0)/2);
			if(i%2!=j%2){
				desl+= 2;
			}

			fprintf(arq, "%e %e %d\n", x, y, phi[i*Nj+j]);
		}
	}
  fclose(arq);
}*/

// Rede quadrada			

			x= i;
			y= j;
		}
	}
	fprintf(arq, "%e %e %d\n", x, y, phi[j*Nj+i]);
	fclose(arq);
}


// Rede triangular			
		
// 			x= 1.0*i+0.5*(j%2);
//			y= 1.0*j*(sqrt(3.0)/2);


void ic(int *phi){ 
	int i, j, cont_inf, trig1;

	cont_inf= 1;
	trig1= 0;
	for(i= 0; i< Ni; i++){
		for(j= 0; j< Nj; j++){
            if(i > Ni/2 && j > Nj/2 && trig1 != 1){
                phi[i*Nj+j]= 2;
				trig1= 1;
            }
            else{
                phi[i*Nj+j]= 1;
            }
		}
	}
}

/*
void recur(int j, int i, int v) {
	if(i >= 0 && j >= 0 && i < L && j < L && phi[j*L+i] == 2){ 
		phi[j*L+i] = v;
		recur((j+1), i, v);
		recur((j-1), i, v);
		recur(j, (i+1), v);
		recur(j, (i-1), v);
	}
}

int count_clusters(void){
	int i, j, cls;
	cls= L*L; 
	for(j= 0; j< L; j++){
		for(i= 0; i< L; i++){
			if(phi[j*L+i] != 2) 
				continue ;
			recur(j, i, ++cls);
		}
	}
	return cls;
}
*/


void vz(int z, int i, int j, int *nb){

	int Nx= Ni;
	int Ny= Nj;
  switch(z){
    case(3):
		nb[0]= ((j-1+Nx)%Nx)*Nx +i;     //top
    nb[1]= ((j+1)%Nx)*Nx    +i;	    //bottom
    if(i%2==j%2){
      nb[2]= j*Nx + (i+1)%Nx;       //horiz-right 
    }else{
      nb[2]= j*Nx + (i-1+Nx)%Nx;    //horiz-left
    }
    break;
    case(4):
      nb[0]= j*Nx+((i+1)%Nx);
      nb[1]= j*Nx+((i-1+Nx)%Nx);
      nb[2]= ((j+1)%Ny)*Nx+i;
      nb[3]= ((j-1+Ny)%Ny)*Nx+i;
    break;
    case(6):
      nb[0]= ((j-1+Nx)%Nx)*Nx +((i+(j%2))%Nx);        //top-right
      nb[1]= j*Nx             +((i+1)%Nx);            //right
      nb[2]= ((j+1)%Nx)*Nx    +((i+(j%2))%Nx);    	  //bot-right
      nb[3]= ((j+1)%Nx)*Nx    +((i-1+Nx+(j%2))%Nx); 	//bot-left
      nb[4]= j*Nx             +((i-1+Nx)%Nx);         //left
      nb[5]= ((j-1+Nx)%Nx)*Nx +((i-1+Nx+(j%2))%Nx);   //top-left
    break;
    case(8):
      nb[0]= ((j-1+Nx)%Nx)*Nx +     i;                //top
      nb[1]= ((j-1+Nx)%Nx)*Nx +((i+1)%Nx);            //top-right
      nb[2]= j*Nx             +((i+1)%Nx);            //right
      nb[3]= ((j+1)%Nx)*Nx    +((i+1)%Nx);          	//bot-right
      nb[4]= ((j+1)%Nx)*Nx    +     i;              	//bottom
      nb[5]= ((j+1)%Nx)*Nx    +((i-1+Nx)%Nx);       	//bot-left
      nb[6]= j*Nx             +((i-1+Nx)%Nx);         //left
      nb[7]= ((j-1+Nx)%Nx)*Nx +((i-1+Nx)%Nx);         //top-left
    break;
    default:
    fprintf(stderr, "Invalid neighbor, check 'z' value\n");
    exit(1);  //kill program with error code 1
  }
}




int main(int argc, char **argv){
	int i, j, x; 
	int ativo, passivo, vizinho, cont_inf, cont_sus, cont_s, cont_i, cont_r, cont_v, trigger= 0;
	int max, t=0;
	int a_ih, a_fh, a_iv, a_fv, key;
	int *nb;		//neighboors
	int z=4;
	char name[100];
	char nome[100];
	double acao;
	FILE *arq;
	
	nb=  (int *) calloc((z), sizeof(int));
		
	if(argc==2){
		gsl_rng_default_seed= atoi(argv[1]); 
	}else{
		gsl_rng_default_seed= time(NULL);
	}
	gsl_rng *w= gsl_rng_alloc(gsl_rng_taus);


	ic(phi);
	cont_sus= 1; // trigger 1
	cont_inf= 1; // trigger 2
	//op(0);
	//
	//
	sprintf(name, "data-%ld.dat", gsl_rng_default_seed);
	arq= fopen(name, "w");
	fclose(arq);
	//while(cont_sus > 0){
	//	while(cont_inf > 0){
			while(t < NG){
				cont_s = cont_i = cont_r = 0;
				for(i= 0; i< L*L; i++){ // verificação do número de indivíduos em cada classe para escrever os dados na pasta plot_curve
					if(phi[i] == 1){
						cont_s++;
					}
					else if(phi[i] == 2){
						cont_i++;
					}
					else if(phi[i] == 3){
						cont_r++;
					}
				}
				sprintf(name, "data-%ld.dat", gsl_rng_default_seed);
				arq= fopen(name, "a");
				fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
				fclose(arq);

				for(x= 0; x< Ni*Nj; x++){
					i= gsl_rng_uniform(w)*Ni;
					j= gsl_rng_uniform(w)*Nj;
					ativo= j*Nj+i;
					
					vz(z, i, j, nb);	
					passivo= nb[(int)(gsl_rng_uniform(w)*z)];
					
					acao= gsl_rng_uniform(w);
					if(acao < beta){
						//infecção
						if(phi[passivo] == 1 && phi[ativo] == 2){
							phi[passivo]= phi[ativo];
						}
					}
					else{
						//recuperação
						if(phi[ativo] == 2){
							phi[ativo]= 3;
						}
					}
				}
				/*cont_sus= cont_inf= 0;
				for(i= 0; i < Ni; i++){
					for(j= 0; j < Nj; j++){
						if(phi[i*Nj+j]==1){
							cont_sus++;
						}
						else if(phi[i*Nj+j] == 2){
							cont_inf++;
						}
					}
				}/*
				if(cont_inf == 0){ // verificação de término prematuro
					if(t>30){ // se a simulação rodou certinho
						cont_sus= 0;
						cont_s = cont_i = cont_r = 0;
						for(i= 0; i< L*L; i++){ // verificação do número de indivíduos em cada classe para escrever os dados na pasta plot_curve
							if(phi[i] == 1){
								cont_s++;
							}
							else if(phi[i] == 2){
								cont_i++;
							}
							else if(phi[i] == 3){
								cont_r++;
							}
						}
						t++;
						sprintf(name, "curve/data-%ld.dat", gsl_rng_default_seed);
						arq= fopen(name, "a");	
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
						//
						//
						//
						//op(t);
						t= NG+1;
						trigger= 1;
					}
					else{ // se a simulação se encerrou prematuramente
						ic(phi);
						t= 0;
						sprintf(name, "curve/data-%ld.dat", gsl_rng_default_seed);
						arq= fopen(name, "w");
						fprintf(arq, "%d %d %d %d\n", t, cont_s, cont_i, cont_r);
						fclose(arq);
					}
				}
				*/
				if(trigger == 0){	
					sprintf(nome, "plot_dat/sir-%d.dat", t);
					arq= fopen(nome, "w");
					for(i= 0; i< Ni; i++){
						for(j= 0; j< Nj; j++){
							fprintf(arq, "%d ", phi[j*Nj+i]);
						}
						fprintf(arq, "\n");
					}
					fclose(arq);
					t++;
					//
					//
					//
				}
			}
		//}
	//}

	
//				Percolação
	/*
	key= 0;
	max= count_clusters();
	for(i= L*L+1; i<= max; i++){ 
		a_ih= a_fh= 0;
		for(j= 0; j< L; j++){
			if(phi[j] == i){ 
				a_ih++; 
			}
			if(phi[(L-1)*L+j] == i){ 
				a_fh++; 
			}
		}

		if(a_ih > 0 && a_fh > 0){
			key= 1;
		}
	}
	
	sprintf(name, "dat/%d_%d/percolacao_%d.dat", repeticoes, L, pr_int);
	if(!(file= fopen(name, "a"))){
		printf("cannot open file teste\n");
		exit(0);
	}
	if(key == 1){
		fprintf(file, "1\n");
		printf("percolou\n");
	}
	else{
		fprintf(file, "0\n");
		printf("não percolou\n");
	}
	fclose(file);
	*/
	gsl_rng_free(w);
	return 0;
}
