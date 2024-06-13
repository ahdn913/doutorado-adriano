#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_statistics.h>

int main(int argc, char **argv){
	if(argc != 2){
		printf("%s input\n", argv[0]);
		exit(0);
	}
	int i, f;
	double x, max;
	double *inp;
	FILE *file;
	char name[100];
	
	if((inp= (double *) calloc(500, sizeof(double))) == NULL){
		printf("cannot allocate inp\n");
		exit(1);
	}
	sprintf(name, "%s", argv[1]);
	if(!(file= fopen(name, "r"))){
		printf("cannot open file %s\n", argv[1]);
		exit(1);
	}
	for(i= 0; i< 500; i++){
		if((fscanf(file, "%lf %lf", &inp[i], &x) != 2)){
			printf("cannot read file %s\n", argv[1]);
			exit(1);
		}
	}
	fclose(file);

	max= 0.0;
	for(i= 1; i < 500; i++){
		if(inp[i] > max){
			max= inp[i];
			f= i-1;
		}
	}

	printf("%d %e\n", f, max);

	free(inp);
	return 0;
}
