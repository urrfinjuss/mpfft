#include "mpfft_header.h"
#ifndef MPFFT_PTHREAD_H
#define MPFFT_PTHREAD_H
#ifdef MPFFT_SERIAL_H
#error "Cannot include both mpfft_serial.h and mpfft_pthread.h in one program"
#endif
#include <pthread.h>
typedef struct {
	unsigned nbits;
	int	nthreads;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpc_ptr	in, out;
	mpc_ptr	W;
} mpfft_plan;

typedef struct {
	mpfr_prec_t precision;
	int tid;
	int stride;
	mpfr_t sum;
	mpfr_t dx;
	mpc_ptr f;
} thread_data;
#endif

//mpfft_pthread.c
extern void mpfft_pthread_init(mpfr_prec_t precision);
extern void mpfft_pthread_destroy_plan(mpfft_plan *plan);
extern void mpfft_pthread_execute(mpfft_plan plan);
extern mpfft_plan mpfft_pthread_create_plan_1d(mpc_ptr out, mpc_ptr in, int nthreads, unsigned nbits, mpfr_prec_t precision, int isign);
extern void mpfft_pthread_example(const int nthreads, mpfr_prec_t precision);

