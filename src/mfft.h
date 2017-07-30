#define MFFT_MAJOR_VERSION	2
#define MFFT_MINOR_VERSION	5
#define mode MPFR_RNDN
#define FFT_FORWARD		-1
#define FFT_BACKWARD	+1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <gmp.h>
#include <mpfr.h>

typedef struct mp_complex mpc, *mpc_ptr;

typedef struct {
	mpfr_prec_t precision;
	int tid;
	int stride;
	mpfr_t sum;
	mpfr_t dx;
	mpc_ptr f;
} thread_data;

struct mp_complex {
	mpfr_t	re;
	mpfr_t	im;
};

typedef struct {
	unsigned nbits;
	int	nthreads;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpc_ptr	in, out;
	mpc_ptr	W;
} mfft_plan;

typedef struct {
	unsigned nbits;
	int	nthreads;
	int	dir;
	double	re, im;
	double complex	*in;
	double complex	*out;
	double complex	*buf;
	double complex	*W;
} base_plan;

// version.c
extern void mfft_version();

//mfft_pthread.c
extern void mfft_pthread_init(mpfr_prec_t precision);
extern void mfft_pthread_destroy_plan(mfft_plan *plan);
extern void mfft_pthread_execute(mfft_plan plan);
extern mfft_plan mfft_pthread_create_plan_1d(mpc_ptr out, mpc_ptr in, int nthreads, unsigned nbits, mpfr_prec_t precision, int isign);
extern void mfft_pthread_example(const int nthreads, mpfr_prec_t precision);

// mfft_serial.c
extern void mfft_init(mpfr_prec_t precision); 
extern void mfft_execute(mfft_plan plan);
extern void mfft_destroy_plan(mfft_plan *plan);
extern mfft_plan mfft_create_plan_1d(mpc_ptr out, mpc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign);

// base_algorithm.c

// init.c
extern void mfft_init(mpfr_prec_t precision);
extern mpc_ptr init_mpc(mpfr_prec_t precision);
extern mpc_ptr init_mpc_array(unsigned N, mpfr_prec_t precision);
extern void mpc_clear(mpc_ptr in);
extern void mpc_clear_array(mpc_ptr in, unsigned N);

// set.c 
extern void mpc_set(mpc_ptr in, mpfr_t x, mpfr_t y);
extern void mpc_set_c(mpc_ptr in, mpc_ptr z);
extern void mpc_set_d(mpc_ptr in, double x, double y);

// bitrev.c
extern void init_bit_operations(mpfr_prec_t precision);
extern void mpfr_bit_reverse(mpc_ptr in, unsigned nbits);
extern void mpfr_bit_reverse_copy(mpc_ptr out, mpc_ptr in, unsigned nbits);
