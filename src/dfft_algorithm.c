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
	if (FFT_TYPE == 1) { // Recursive
		double complex tmp;
		plan->tmp = malloc((1<<nbits)*sizeof(double complex));
		for (unsigned s = 0; s < nbits; s++) {
			tmp = M_PI/(1<<s);
			plan->W[s] = cexp(1.I*tmp*isign);
 		}	
		return *plan;
	} else if (FFT_TYPE == 2) { // Danielson-Lanczos
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


//fft_recursive(plan.out, plan.tmp, n, 1);
//static double complex t  = 1;
//void fft_recursive(double complex *out, double complex *tmp, int n, int stride) {
static int S = 0;
void fft_recursive(double complex *out, double complex *tmp, double complex *tw, int n, int stride) {
	if (stride < n) {
		fft_recursive(tmp, out, tw, n, 2*stride);
		fft_recursive(tmp + stride, out + stride, tw, n, 2*stride);
		double complex t0 = cexp(-2.I*M_PI*stride/n);
		double complex t = 1;
		//int ss = 0;
		for (int j = 0; j < n; j += 2*stride) {
			//double complex t = cexp(-1.I*M_PI*j/n)*out[j+stride];
			//double complex t = cexp(-1.I*M_PI*j/n);
			out[j/2]     = tmp[j] + t*tmp[j+stride];
			out[(j+n)/2] = tmp[j] - t*tmp[j+stride];
			//t = t*tw[S];
			t = t*t0;
		}
		//printf("ss = %d\n", S);
		//S++;
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
		memcpy(plan.tmp, plan.in, n*sizeof(double complex));
		memcpy(plan.out, plan.in, n*sizeof(double complex));
		S = 0;
		fft_recursive(plan.out, plan.tmp, plan.W, n, 1);
		//fft_recursive(plan.out, plan.tmp, n, 1);
	} else if (plan.FFT_TYPE == 2) {
		fft_danielson_lanczos(plan);
	}
}

void fft_destroy_plan(fft_plan plan) {
	free(plan.tmp);
}

