#ifndef MFFT_HEADER_H
#include "mfft_header.h"
#endif

// mfft_serial.c
extern void mfft_init(mpfr_prec_t precision); 
extern void mfft_execute(mfft_plan plan);
extern void mfft_destroy_plan(mfft_plan *plan);
extern mfft_plan mfft_create_plan_1d(mpc_ptr out, mpc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign);

