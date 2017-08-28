#include "mpfft_header.h"
#ifndef MPFFT_SERIAL_H
#define MPFFT_SERIAL_H
#ifdef MPFFT_PTHREAD_H
#error "Cannot include both mpfft_serial.h and mpfft_pthread.h in one program"
#endif
typedef struct {
	unsigned nbits;
	unsigned int inplace;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpfc_ptr	in, out;
	mpfc_ptr	W;
} mpfft_plan;
#endif

// mpfft_serial.c
extern void mpfft_init(mpfr_prec_t precision); 
extern void mpfft_execute(mpfft_plan plan);
extern void mpfft_destroy_plan(mpfft_plan plan);
extern mpfft_plan mpfft_create_plan_1d(mpfc_ptr out, mpfc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign);

