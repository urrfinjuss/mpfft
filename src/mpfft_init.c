#include "mpfft_header.h"

mpfc_ptr init_mpfc(mpfr_prec_t precision) {
	mpfc_ptr in = malloc(sizeof(mpfc_t));
	mpfr_init2(in->re, precision);
	mpfr_init2(in->im, precision);
	mpfr_set_ui(in->re, 0, MODE);
	mpfr_set_ui(in->im, 0, MODE);
	return in;
}

mpfc_ptr init_mpfc_array(unsigned N, mpfr_prec_t precision) {
	mpfc_ptr in = malloc(N*sizeof(mpfc_t));
	for (unsigned j = 0; j < N; j++) {
		mpfr_init2(in[j].re, precision);
		mpfr_init2(in[j].im, precision);
		mpfr_set_ui(in[j].re, 0, MODE);
		mpfr_set_ui(in[j].im, 0, MODE);
	}
	return in;
}

void mpfc_clear(mpfc_ptr in) {
	mpfr_clear(in->re);
	mpfr_clear(in->im);
}

void mpfc_clear_array(mpfc_ptr in, unsigned N) {
	for (unsigned j = 0; j < N; j++) {
		mpfr_clear(in[j].re);
		mpfr_clear(in[j].im);
	}
}

