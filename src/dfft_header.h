#include "header.h"
#include <math.h>
#include <complex.h>

typedef struct {
	int isign;
	int FFT_TYPE;
	unsigned nbits;
	double complex *tmp;
	double complex *in;
	double complex *out;
	double complex *W;
} fft_plan;


// dfft_algorithm.c
extern fft_plan fft_create_plan_1d(double complex *in, double complex *out, unsigned nbits, int isign, int FFT_TYPE);
extern void fft_recursive(double complex *buf, double complex *out, int n, int stride);
extern void fft_execute(fft_plan plan);
extern void fft_destroy_plan(fft_plan plan);

// dfft_bitrev.c
void dfft_bit_reverse(double complex *in, unsigned nbits);
void dfft_bit_reverse_copy(double complex *out, double complex *in, unsigned nbits);
