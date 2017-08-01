#include "dfft_header.h"

static double complex temp;

void dfft_bit_reverse(double complex *in, unsigned nbits) {
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
			temp = in[n];
			in[n] = in[nreversed];
			in[nreversed] = temp;
		}
	}
}

void dfft_bit_reverse_copy(double complex *out, double complex *in, unsigned nbits){
	if (out == in) printf("Error: In-Place Transform!\n");
	for (int j = 0; j < 1<<nbits; j++) {
		out[j] = in[j];
	}
	dfft_bit_reverse(out, nbits);
}
