#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define Nx 100
#define Ny Nx
#define size_x 1000						//Size x of image
#define size_y 1000						//Size y of image
#define max_x 99.5						//Max value of data in x
#define max_y	85.73651				//Max value of data in y

#define N_SITES Nx*Ny
#define	mult_x 1.0*size_x/max_x			//size_x*Max_x				
#define	mult_y 1.0*size_y/max_y 		//size_x*Max_x

double site[N_SITES][2];
unsigned char rgb[N_SITES][3];

//int size_x = 1000, size_y = 1000;

inline double sq2(double x, double y)
{
	return x * x + y * y;
}

#define for_k for (k = 0; k < N_SITES; k++)
int nearest_site(double x, double y)
{
	int k, ret = 0;
	double d, dist = 0;
	for_k {
		d = sq2(x - site[k][0], y - site[k][1]);
		if (!k || d < dist) {
			dist = d, ret = k;
		}
	}
	return ret;
}

/* see if a pixel is different from any neighboring ones */
int at_edge(int *color, int y, int x)
{
	int i, j, c = color[y * size_x + x];
	for (i = y - 1; i <= y + 1; i++) {
		if (i < 0 || i >= size_y) continue;

		for (j = x - 1; j <= x + 1; j++) {
			if (j < 0 || j >= size_x) continue;
			if (color[i * size_x + j] != c) return 1;
		}
	}
	return 0;
}

#define AA_RES 4 /* average over 4x4 supersampling grid */
void aa_color(unsigned char *pix, int y, int x)
{
	int i, j, n;
	double r = 0, g = 0, b = 0, xx, yy;
	for (i = 0; i < AA_RES; i++) {
		yy = y + 1. / AA_RES * i + .5;
		for (j = 0; j < AA_RES; j++) {
			xx = x + 1. / AA_RES * j + .5;
			n = nearest_site(xx, yy);
			r += rgb[n][0];
			g += rgb[n][1];
			b += rgb[n][2];
		}
	}
	pix[0] = r / (AA_RES * AA_RES);
	pix[1] = g / (AA_RES * AA_RES);
	pix[2] = b / (AA_RES * AA_RES);
}

#define for_i for (i = 0; i < size_y; i++)
#define for_j for (j = 0; j < size_x; j++)
void gen_map()
{
	int i, j, k;
	int *nearest = malloc(sizeof(int) * size_y * size_x);
//	unsigned char *ptr, *buf, color;
	unsigned char *ptr, *buf;

	ptr = buf = malloc(3 * size_x * size_y);
	for_i for_j nearest[i * size_x + j] = nearest_site(j, i);

	for_i for_j {
		if (!at_edge(nearest, i, j))
			memcpy(ptr, rgb[nearest[i * size_x + j]], 3);
		else	/* at edge, do anti-alias rastering */
			aa_color(ptr, i, j);
		ptr += 3;
	}

	/* draw sites */
	for (k = 0; k < N_SITES; k++) {
//		color = (rgb[k][0]*.25 + rgb[k][1]*.6 + rgb[k][2]*.15 > 80) ? 0 : 255;

		for (i = site[k][1] - 1; i <= site[k][1] + 1; i++) {
			if (i < 0 || i >= size_y) continue;

			for (j = site[k][0] - 1; j <= site[k][0] + 1; j++) {
				if (j < 0 || j >= size_x) continue;

//				ptr = buf + 3 * (i * size_x + j);
//				ptr[0] = ptr[1] = ptr[2] = color;
			}
		}
	}

	printf("P6\n%d %d\n255\n", size_x, size_y);
	fflush(stdout);
	fwrite(buf, size_y * size_x * 3, 1, stdout);
}

int main(int argc, char *argv[])
{
	if(argc != 2){
		printf("Uso: %s <sir-X.dat>\n", argv[0]);	
		printf("<sir-X.dat>: Valor de X em sir-X.dat\n");	
		return 1;		
	}

	int k; 
	double data[N_SITES][3];
	char arq[50];

	sprintf(arq, "sir-%d.dat", atoi(argv[1]));
	FILE *arqv= fopen(arq, "r");

	if(arqv == NULL){
		printf("Falha ao abrir o arquivo %s\n", arq);
		return 1;
	}

	for_k {
		if(fscanf(arqv, "%lf %lf %lf", &data[k][0], &data[k][1], &data[k][2]) != 3){
			printf("Falha ao abrir o arquivo %s\n", arq);
			return 1;
		}
//		printf("-------------------\n");		
//		printf("data[%d][1] = %lf\n", k, data[k][0]);		
//		printf("data[%d][2] = %lf\n", k, data[k][1]);		
//		printf("data[%d][3] = %lf\n", k, data[k][2]);	
	}

	//verificação Escalonamento
	if(size_x < Nx || size_y < Ny){
		printf("Erro: Size x/y menores que Nx/Ny");
		return 1;		
	}

	for_k {
		site[k][0] = 				mult_x*data[k][0];
		site[k][1] = size_y-mult_y*data[k][1];

		switch ((int)data[k][2]){
			case 1:
				rgb[k][0] = 80;
				rgb[k][1] = 255;
				rgb[k][2] = 80;
				break;
	
			case 2:
				rgb[k][0] = 255;
				rgb[k][1] = 80;
				rgb[k][2] = 80;
				break;
			
			default:
				rgb[k][0] = 80;
				rgb[k][1] = 80;
				rgb[k][2] = 255;
				break;
		}
		
//		printf("------------------\n");
//		printf("site [k][0]= %d \n", (int)site[k][0]);
//		printf("site [k][1]= %d \n", (int)site[k][1]);
//		printf("rgb  [k][0]= %d \n",  rgb[k][0]);
//		printf("rgb  [k][1]= %d \n",  rgb[k][1]);
//		printf("rgb  [k][2]= %d \n",  rgb[k][2]);
	}

	gen_map();
	return 0;
}
