//#ifndef MFFT_HEADER_H
//#include "mfft_header.h"
//#endif
#include "base_fft.h"


fft_plan fft_create_plan(double complex *in, double complex *out, int n, int isign) {
	fft_plan *plan = malloc(sizeof(fft_plan));
	plan->n = n;
	plan->isign = isign;
	plan->in = in;
	plan->out = out;
	plan->tmp = malloc(n*sizeof(double complex));
	memcpy(plan->tmp, plan->in, plan->n*sizeof(double complex));
	memcpy(plan->out, plan->in, plan->n*sizeof(double complex));
	return *plan;
}


void fft_recursive(double complex *buf, double complex *out, int n, int stride) {
	if (stride < n) {
		fft_recursive(out, buf, n, 2*stride);
		fft_recursive(out + stride, buf + stride, n, 2*stride);
		for (int j = 0; j < n; j += 2*stride) {
			double complex t = cexp(-1.I*M_PI*j/n)*out[j+stride];
			buf[j/2]     = out[j] + t;
			buf[(j+n)/2] = out[j] - t;
		}
	}
}

void call_fft(fft_plan plan) {
	fft_recursive(plan.out, plan.tmp, plan.n, 1);
}

