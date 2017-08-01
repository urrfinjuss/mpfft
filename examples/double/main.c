#include <stdio.h>
#include <string.h>
#include <dfft_header.h>
#include <math.h>
#include <complex.h>
#include <time.h>

#define DFFT_RECURSIVE 1
#define DFFT_DANIELSON_LANCZOS 2

int main() {
	mpfft_version();
	double complex *y, *x;
	double t;
	unsigned nbits = 10;
	int n = 1<<nbits;
	FILE *fh = fopen("init.txt","w");
	fprintf(fh, "# 1. x 2. y\n\n");

	x = malloc(n*sizeof(double complex));
	y = malloc(n*sizeof(double complex));
	for (int j = 0; j < n ; j++) {
		t = 2.*M_PI*j/n;
		x[j] = sin(t);
		fprintf(fh, "%17.12e\t%17.12e\t%17.12e\n", t, creal(x[j]), cimag(x[j]));
	}
	fclose(fh);
	fft_plan plan = fft_create_plan_1d(x, y, nbits, 1, DFFT_DANIELSON_LANCZOS);
	fft_execute(plan);

	int k;
	fh = fopen("fft.txt","w");
	fprintf(fh, "# 1. k 2. re z_k 3. im z_k\n\n");
	for (int j = 0; j < n; j++) {
		k = j;
		if (j > n/2) k = j - n;
		fprintf(fh, "%5d\t%17.12e\t%17.12e\n", k, creal(y[j]), cimag(y[j]));
	}
	fclose(fh);


	return 0;
}
