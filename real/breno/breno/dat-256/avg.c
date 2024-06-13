#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics.h>

int main(int argc, char **argv){
	if(argc != 5){
		printf("\nplease, type:\n\n %s name_of_file n_columns n_rows n_files\n\n", argv[0]);
		exit(0);
	}
	int i, j, k, NC, NR, NF;
	double *avg, *inp;
	FILE *file;
	char name[100];
	
	NC= atoi(argv[2]);
	NR= atoi(argv[3]);
	NF= atoi(argv[4]);

	if((avg= (double *) calloc(NF, sizeof(double))) == NULL){
		printf("cannot allocate avg\n");
		exit(1);
	}
	if((inp= (double *) calloc(NR*NC*NF, sizeof(double))) == NULL){
		printf("cannot allocate inp\n");
		exit(1);
	}
	for(k= 0; k< NF; k++){
		sprintf(name, "%s-%d.dat", argv[1], k+1);
		if(!(file= fopen(name, "r"))){
			printf("cannot open file %s-%d.dat\n", argv[1], k+1);
			exit(1);
		}
		for(j= 0; j< NR; j++){
			for(i= 0; i< NC; i++){
				if((fscanf(file, "%lf", &inp[(k*NR+j)*NC+i]) != 1)){
					printf("cannot read file %s-%d.dat\n", argv[1], k+1);
					exit(1);
				}
			}
		}
		fclose(file);
	}

	sprintf(name, "%s-avg.dat", argv[1]);
	if(!(file= fopen(name, "w"))){
		printf("cannot open file %s-avg.dat\n", argv[1]);
		exit(1);
	}
	for(j= 0; j< NR; j++){
		for(i= 0; i< NC; i++){
			for(k= 0; k< NF; k++){
				avg[k]= inp[(k*NR+j)*NC+i];
			}
			fprintf(file, "%+e ", gsl_stats_mean(avg, 1, NF));
		}
		fprintf(file, "\n");
	}
	fclose(file);

	free(inp);
	free(avg);
	return 0;
}
