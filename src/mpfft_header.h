#ifndef MPFFT_HEADER_H
#define MODE MPFR_RNDN
#define MPFFT_HEADER_H
#include "header.h"
#include <gmp.h>
#include <mpfr.h>

typedef struct mp_complex {
	mpfr_t	re;
	mpfr_t	im;
} mpfc_t, *mpfc_ptr;

#endif

// mpfr_version.c
extern void mpfr_my_version();

// init.c
extern void mpfft_init(mpfr_prec_t precision);
extern mpfc_ptr init_mpfc(mpfr_prec_t precision);
extern mpfc_ptr init_mpfc_array(unsigned N, mpfr_prec_t precision);
extern void mpfc_clear(mpfc_ptr in);
extern void mpfc_clear_array(mpfc_ptr in, unsigned N);

// set.c 
extern void mpfc_set(mpfc_ptr in, mpfr_t x, mpfr_t y);
extern void mpfc_set_c(mpfc_ptr in, mpfc_ptr z);
extern void mpfc_set_d(mpfc_ptr in, double x, double y);

// bitrev.c
extern void init_bit_operations(mpfr_prec_t precision);
extern void mpfr_bit_reverse(mpfc_ptr in, unsigned nbits);
extern void mpfr_bit_reverse_copy(mpfc_ptr out, mpfc_ptr in, unsigned nbits);

