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
	unsigned nbits = 16;
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
	
	int ntimes = 2;
	struct timespec begin, end;
	double elapsed;

	
	fft_plan plan = fft_create_plan_1d(x, y, nbits, 1, DFFT_DANIELSON_LANCZOS);
	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (int j = 0; j < ntimes; j++) {
		fft_execute(plan);
	}	
	clock_gettime(CLOCK_MONOTONIC, &end);
	elapsed = end.tv_sec - begin.tv_sec;
	elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
	printf("Danielson-Lanczos FFT algorithm:\t");
	printf("avg time for size %d FFT %.2f ms\n", n, 1e+3*elapsed/ntimes);	
	
	
	fft_plan plan1 = fft_create_plan_1d(x, y, nbits, 1, DFFT_RECURSIVE);
	clock_gettime(CLOCK_MONOTONIC, &begin);
	for (int j = 0; j < ntimes; j++) {
		fft_execute(plan1);
	}	
	clock_gettime(CLOCK_MONOTONIC, &end);
	elapsed = end.tv_sec - begin.tv_sec;
	elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
	printf("Recursive FFT algorithm:\t\t");
	printf("avg time for size %d FFT %.2f ms\n", n, 1e+3*elapsed/ntimes);	
	

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
