#define MFFT_MAJOR_VERSION	2
#define MFFT_MINOR_VERSION	2
#define mode MPFR_RNDN
#define FFT_FORWARD		-1
#define FFT_BACKWARD	+1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gmp.h>
#include <mpfr.h>

typedef struct mp_complex mpc, *mpc_ptr;
typedef struct mfft_naive naive_plan;
typedef struct fft_plan mfft_plan;

struct mp_complex {
	mpfr_t	re;
	mpfr_t	im;
};

struct fft_plan {
	unsigned nbits;
	int	dir;
	mpfr_prec_t prec;
	mpfr_t	re, im;
	mpc_ptr	in, out;
	mpc_ptr	W;
};

struct mfft_naive {
	unsigned int N;
	mpc_ptr	w;
	mpc_ptr	in;
	mpc_ptr	out;
	mpfr_t	tmp1;
	mpfr_t	tmp2;
};

// version.c
extern void mfft_version();

// mfft.c
extern void init_mfft(mpfr_prec_t precision); 
extern mfft_plan mfft_create_plan_1d(mpc_ptr out, mpc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign);
extern void mfft_execute(mfft_plan plan);
extern void mfft_destroy_plan(mfft_plan *plan);

// dfft.c
extern naive_plan init_mfft_naive_plan(mpc_ptr out, mpc_ptr in, unsigned N, mpfr_prec_t precision, int isign);
extern void mfft_naive_execute(naive_plan p);

// init.c
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
