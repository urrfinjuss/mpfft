#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <dfft_header.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>
#include <limits.h>

#define NUM_THREADS 4

#define DFFT_RECURSIVE 1
#define DFFT_DANIELSON_LANCZOS 2
#define DFFT_FFTW	3
#define DDFT_NAIVE	4

static double complex *x, *y;

unsigned int intlog2 (int val) {
    if (val == 0) return UINT_MAX;
    if (val == 1) return 0;
    unsigned int ret = 0;
    while (val > 1) {
        val >>= 1;
        ret++;
    }
    return ret;
}

void call_init_test_arrays(int N) {
	FILE *fh = fopen("init.txt","w");
	x = fftw_malloc(N*sizeof(double complex));
	y = fftw_malloc(N*sizeof(double complex));
	memset(y, 0, N*sizeof(double complex));
	fprintf(fh, "# 1. x 2. y\n# N = %d\n\n", N);
	for (int j = 0; j < N; j++) {
		x[j] = sin(2.*M_PI*j/N);
		fprintf(fh, "%19.12e\t%19.12e\t%19.12e\n", 2.*M_PI*j/N, creal(x[j]), cimag(x[j]));
	}
	fclose(fh);
}

void call_check_result(int N) {
	char fname[80];
	int k = 0;
	sprintf(fname, "fft_n%d.txt", N);
	FILE *fh = fopen(fname, "w");
	fprintf(fh, "# 1. k 2. re z_k 3. im z_k\n# N = %d\n\n", N);
	for (int j = 0; j < N; j++) {
		k = j;
		if (j > N/2) k = j - N;
		fprintf(fh, "%5d\t%19.12e\t%19.12e\n", k, creal(y[j]), cimag(y[j]));
	}
	fclose(fh);
}

int main(int argc, char **argv) {
	mpfft_version();
	if (argc != 3) {
		printf("Usage: %s fft_type size\n", argv[0]);
		printf("Available FFT types: recursive, danielson_lanczos, fftw, nsquared\n");
		exit(0);
	}
	int N = atoi(argv[2]);
	int nbits = intlog2(N);
	struct timespec begin, end;
	double elapsed;
	int ntimes = 1024;
	call_init_test_arrays(N);
	if (!strcmp(argv[1], "recursive")) {
		printf("Using recursive algorithm N = %6d: ", N);
		fft_plan plan = fft_create_plan_1d(x, y, nbits, 1, DFFT_RECURSIVE);
		clock_gettime(CLOCK_MONOTONIC, &begin);
		for (int j = 0; j < ntimes; j++) fft_execute(plan);
		clock_gettime(CLOCK_MONOTONIC, &end);
		elapsed = end.tv_sec - begin.tv_sec;
		elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
		fft_destroy_plan(plan);
		call_check_result(N);	
		printf("\t%.2f ms  ( %8.2f ns)\n", 1e+3*elapsed/ntimes, 1e+6*elapsed/ntimes);
		exit(0);
	} else if (!strcmp(argv[1],"danielson_lanczos")) {
		printf("Using Danielson-Lanczos algorithm N = %6d: ", N);
		fft_plan plan = fft_create_plan_1d(x, y, nbits, 1, DFFT_DANIELSON_LANCZOS);
		clock_gettime(CLOCK_MONOTONIC, &begin);
		for (int j = 0; j < ntimes; j++) fft_execute(plan);
		clock_gettime(CLOCK_MONOTONIC, &end);
		elapsed = end.tv_sec - begin.tv_sec;
		elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
		fft_destroy_plan(plan);
		call_check_result(N);	
		printf("\t%.2f ms  ( %8.2f ns)\n", 1e+3*elapsed/ntimes, 1e+6*elapsed/ntimes);
		exit(0);
	} else if (!strcmp(argv[1], "fftw")) {
		printf("Using FFTW with FFTW_ESTIMATE flag N = %6d: ", N);
		fftw_plan plan = fftw_plan_dft_1d(N, x, y, FFTW_FORWARD, FFTW_ESTIMATE);
		clock_gettime(CLOCK_MONOTONIC, &begin);
		for (int j = 0; j < ntimes; j++) fftw_execute(plan);
		clock_gettime(CLOCK_MONOTONIC, &end);
		elapsed = end.tv_sec - begin.tv_sec;
		elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
		fftw_destroy_plan(plan);
		call_check_result(N);	
		printf("\t%.2f ms  ( %8.2f ns)\n", 1e+3*elapsed/ntimes, 1e+6*elapsed/ntimes);
		exit(0);
        } else if (!strcmp(argv[1], "nsquared")) {
		printf("Running naive algorithm N^2 in serial N = %6d: ", N);
		fft_plan plan = dft_create_plan_1d(x, y, nbits, 1, NUM_THREADS, DDFT_NAIVE);
		clock_gettime(CLOCK_MONOTONIC, &begin);
		for (int j = 0; j < ntimes; j++) dft_execute(plan);
		clock_gettime(CLOCK_MONOTONIC, &end);
		elapsed = end.tv_sec - begin.tv_sec;
		elapsed += (end.tv_nsec - begin.tv_nsec) / 1000000000.0;
		fft_destroy_plan(plan);
		call_check_result(N);	
		printf("\t%.2f ms  ( %8.2f ns)\n", 1e+3*elapsed/ntimes, 1e+6*elapsed/ntimes);
		exit(0);	
	} else {
		printf("Unknown FFT algorithm.\n");
		exit(0);
	}
	

	return 0;
}
