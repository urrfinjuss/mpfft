#include "mpfft_header.h"

mpc_ptr init_mpc(mpfr_prec_t precision) {
	mpc_ptr in = malloc(sizeof(mpc));
	mpfr_init2(in->re, precision);
	mpfr_init2(in->im, precision);
	return in;
}

mpc_ptr init_mpc_array(unsigned N, mpfr_prec_t precision) {
	mpc_ptr in = malloc(N*sizeof(mpc));
	for (unsigned j = 0; j < N; j++) {
		mpfr_init2(in[j].re, precision);
		mpfr_init2(in[j].im, precision);
	}
	return in;
}

void mpc_clear(mpc_ptr in) {
	mpfr_clear(in->re);
	mpfr_clear(in->im);
}

void mpc_clear_array(mpc_ptr in, unsigned N) {
	for (unsigned j = 0; j < N; j++) {
		mpfr_clear(in[j].re);
		mpfr_clear(in[j].im);
	}
}

