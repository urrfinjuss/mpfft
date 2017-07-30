#include <stdio.h>
#include <string.h>
#include <stdlib.h>
//#include <mfft.h>
#include <math.h>
#include <complex.h>
#include <time.h>

typedef struct {
	int n;
	int isign;
	double complex *tmp;
	double complex *in;
	double complex *out;
} fft_plan;


// base_algorithms.c
extern fft_plan fft_create_plan(double complex *in, double complex *out, int n, int isign);
extern void fft_recursive(double complex *buf, double complex *out, int n, int stride);
extern void call_fft(fft_plan plan);
extern void fft_destroy_plan(fft_plan plan);
