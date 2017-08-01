#ifndef MFFT_HEADER_H
#include "mfft_header.h"
#endif

#ifndef MFFT_H
#define MFFT_H

/*
typedef struct {
	unsigned nbits;
	int	nthreads;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpc_ptr	in, out;
	mpc_ptr	W;
} mfft_plan;
*/

typedef struct {
	mpfr_prec_t precision;
	int tid;
	int stride;
	mpfr_t sum;
	mpfr_t dx;
	mpc_ptr f;
} thread_data;
#endif

//mfft_pthread.c
extern void mfft_pthread_init(mpfr_prec_t precision);
extern void mfft_pthread_destroy_plan(mfft_plan *plan);
extern void mfft_pthread_execute(mfft_plan plan);
extern mfft_plan mfft_pthread_create_plan_1d(mpc_ptr out, mpc_ptr in, int nthreads, unsigned nbits, mpfr_prec_t precision, int isign);
extern void mfft_pthread_example(const int nthreads, mpfr_prec_t precision);

