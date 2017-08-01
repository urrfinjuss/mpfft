#include "mpfft_header.h"
#ifndef MPFFT_SERIAL_H
#define MPFFT_SERIAL_H
#ifdef MPFFT_PTHREAD_H
#error "Cannot include both mpfft_serial.h and mpfft_pthread.h in one program"
#endif
typedef struct {
	unsigned nbits;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpc_ptr	in, out;
	mpc_ptr	W;
} mpfft_plan;
#endif

// mpfft_serial.c
extern void mpfft_init(mpfr_prec_t precision); 
extern void mpfft_execute(mpfft_plan plan);
extern void mpfft_destroy_plan(mpfft_plan plan);
extern mpfft_plan mpfft_create_plan_1d(mpc_ptr out, mpc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign);

