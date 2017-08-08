#include "header.h"
#include <math.h>
#include <complex.h>

typedef struct {
	const int isign;
	const int FFT_TYPE;
	const unsigned nbits;
        unsigned nthreads;
	double complex *tmp;
	double complex *in;
	double complex *out;
	double complex *W;
} fft_plan;


// dfft_algorithm.c
extern fft_plan fft_create_plan_1d(double complex *in, double complex *out, const unsigned nbits, const int isign, const int FFT_TYPE);
extern void fft_recursive(double complex *buf, double complex *out, double complex *tw, const int n, const int stride);
extern void fft_execute(fft_plan plan);
extern void fft_destroy_plan(fft_plan plan);

// dfft_bitrev.c
void dfft_bit_reverse(double complex *in, unsigned nbits);
void dfft_bit_reverse_copy(double complex *out, double complex *in, const unsigned nbits);

// ddft_algorithm.c
extern void dft_execute(fft_plan p);
extern void dft_execute_pthreads(fft_plan p);
fft_plan dft_create_plan_1d(double complex *in, double complex *out, const unsigned nbits, const int isign, const unsigned nthreads, const int DFT_TYPE);
extern void dft_destroy_plan(fft_plan plan);
