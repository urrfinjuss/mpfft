#ifndef MFFT_HEADER_H
#include "mfft_header.h"
#endif

void mpc_set_c(mpc_ptr in, mpc_ptr z){
	mpfr_set(in->re, z->re, MODE);
	mpfr_set(in->im, z->im, MODE);
}

void mpc_set(mpc_ptr in, mpfr_t x, mpfr_t y){
	mpfr_set(in->re, x, MODE);
	mpfr_set(in->im, y, MODE);
}

void mpc_set_d(mpc_ptr in, double x, double y){
	mpfr_set_d(in->re, x, MODE);
	mpfr_set_d(in->im, y, MODE);
}

void bit_reversal(mpc_ptr in, int N, mpfr_prec_t precision) {
	mpfr_t tempr, tempi;
	mpfr_init2(tempr, precision);
	mpfr_init2(tempi, precision);
	int j = 1;
	for (int i = 1; i < N; i++) {
		if (j > i) {
			mpfr_set(tempr, in[j].re, MODE);
			mpfr_set(tempi, in[j].im, MODE);
			mpfr_set(in[j].re, in[i].re, MODE);
			mpfr_set(in[j].im, in[i].im, MODE);
			mpfr_set(in[i].re, tempr, MODE);
			mpfr_set(in[i].im, tempi, MODE);
		}
		int m = N/2;
		while ((m >= 2)&&(j > m)) {
			j = j - m;
			m = m/2;
		}
		j = j + m;
	}
}
