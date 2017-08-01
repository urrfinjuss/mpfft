#ifndef MFFT_HEADER_H
#include "mfft_header.h"
#endif
#include <math.h>
#include <complex.h>

typedef struct {
	int n;
	int isign;
	double complex *tmp;
	double complex *in;
	double complex *out;
} fft_plan;


// base_algorithms.c
extern fft_plan fft_create_plan_1d(double complex *in, double complex *out, int n, int isign);
extern void fft_recursive(double complex *buf, double complex *out, int n, int stride);
extern void fft_execute(fft_plan plan);
extern void fft_destroy_plan(fft_plan plan);
