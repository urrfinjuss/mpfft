#include "dfft_header.h"

static double complex t, u, w;

fft_plan fft_create_plan_1d(double complex *in, double complex *out, unsigned nbits, int isign, int FFT_TYPE) {
	fft_plan *plan = malloc(sizeof(fft_plan));
	plan->in = in;
	plan->out = out;
	plan->nbits = nbits;
	plan->isign = isign;
	plan->W = malloc(nbits*sizeof(double complex));
	plan->FFT_TYPE = FFT_TYPE;
	if (FFT_TYPE == 1) {
		printf("Recursive FFT algorithm.\n");
		plan->tmp = malloc((1<<nbits)*sizeof(double complex));
		memcpy(plan->tmp, plan->in, (1<<nbits)*sizeof(double complex));
		memcpy(plan->out, plan->in, (1<<nbits)*sizeof(double complex));
		return *plan;
	} else if (FFT_TYPE == 2) {
		printf("Danielson-Lanczos FFT algorithm.\n");
		double complex tmp;
		for (unsigned s = 0; s < nbits; s++) {
			tmp = M_PI/(1<<s);
			plan->W[s] = cexp(1.I*tmp*isign);
 		}	
		return *plan;
	} else {
		printf("Unknown FFT flag.\n");
		exit(0);
	}
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

void fft_danielson_lanczos(fft_plan plan) {
	dfft_bit_reverse_copy(plan.out, plan.in, plan.nbits);
	int n = 1<<plan.nbits;
	for (unsigned s = 0; s < plan.nbits; s++) {
		int m = 1 << (s+1);
		for (int k = 0; k < n; k = k+m) {
			w = 1;
			for (int j = 0; j < m/2; j++) {
				t = w*plan.out[k+j+m/2];
				u = plan.out[k+j];
				plan.out[k+j] = u + t;
				plan.out[k+j+m/2] = u - t;
				w = w*plan.W[s];
			}
		}
	}
}

void fft_execute(fft_plan plan) {
	int n = 1 << plan.nbits;
	if (plan.FFT_TYPE == 1) {
		fft_recursive(plan.out, plan.tmp, n, 1);
	} else if (plan.FFT_TYPE == 2) {
		fft_danielson_lanczos(plan);
	}
}

void fft_destroy_plan(fft_plan plan) {
	free(plan.tmp);
}

