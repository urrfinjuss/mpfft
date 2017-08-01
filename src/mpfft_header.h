#ifndef MPFFT_HEADER_H
#define MODE MPFR_RNDN
#define MPFFT_HEADER_H
#include "header.h"
#include <gmp.h>
#include <mpfr.h>

typedef struct mp_complex {
	mpfr_t	re;
	mpfr_t	im;
} mpc, *mpc_ptr;

#endif

// mpfr_version.c
extern void mpfr_my_version();

// init.c
extern void mpfft_init(mpfr_prec_t precision);
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

