#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

int main(int argc, char **argv) {
	int f, t, N;
	double a, b, c;
	FILE *file;

	if(argc != 4) {
		printf("%s input output N\n", argv[0]);
		exit(1);
	}

	N = atoi(argv[3]);

	if(!(file = fopen(argv[1], "r"))){
		printf("Cannot open file %s\n", argv[1]);
		exit(1);
	}
	fftw_complex *g_t = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);
	fftw_complex *g_f = (fftw_complex *) fftw_malloc(sizeof(fftw_complex)*N);
	for(t = 0; t < N; t++) {
		if(fscanf(file, "%lf %lf %lf", &a, &b, &c) != 3){
			printf("Cannot read from file %s\n", argv[1]);
			exit(1);
		}
		g_t[t] = a + 0.0*I;
	}
	fclose(file);

	fftw_plan FTF = fftw_plan_dft_1d(N, g_t, g_f, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(FTF);
	fftw_destroy_plan(FTF);

	g_f[0]/= N;
	for(f = 1; f < N; f++){
		g_f[f]/= (N/2);
	}
	if(!(file = fopen(argv[2], "w"))){
		printf("Cannot open file %s\n", argv[2]);
		exit(1);
	}
	for(f = 0; f < N/2; f++){
		fprintf(file, "%e\n", cabs(g_f[f]*g_f[f]));
	}
	fclose(file);

	fftw_free(g_t);
	fftw_free(g_f);
	fftw_cleanup();

	return 0;
}
