#include "mpfft_header.h"

static mpfr_t tempr, tempi;

void init_bit_operations(mpfr_prec_t precision) {
	mpfr_init2(tempr, precision);
	mpfr_init2(tempi, precision);
}

void mpfr_bit_reverse(mpfc_ptr in, unsigned nbits) {
	for (int n  = 1; n < 1<<nbits; n++) {
		int nreversed = n;
		int count = nbits-1;
		for (int nforward = n>>1; nforward; nforward >>= 1) {
			nreversed <<= 1;
			nreversed |= nforward & 1;
			count--;
		}
		nreversed <<= count;
		nreversed &= (1<<nbits) - 1;
		if (n < nreversed) {
			mpfr_set(tempr, in[n].re, MODE);
			mpfr_set(tempi, in[n].im, MODE);
			mpfr_set(in[n].re, in[nreversed].re, MODE);
			mpfr_set(in[n].im, in[nreversed].im, MODE);
			mpfr_set(in[nreversed].re, tempr, MODE);
			mpfr_set(in[nreversed].im, tempi, MODE);
		}
	}
}

void mpfr_bit_reverse_copy(mpfc_ptr out, mpfc_ptr in, unsigned nbits) {
	if (out == in) printf("Error: In-Place Transform!\n");
	for (int j = 0; j < 1<<nbits; j++) {
		mpfr_set(out[j].re, in[j].re, MODE);
		mpfr_set(out[j].im, in[j].im, MODE);
	}
	mpfr_bit_reverse(out, nbits);
}
