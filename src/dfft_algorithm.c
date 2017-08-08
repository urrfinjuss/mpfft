#include "dfft_header.h"

static double complex t, w;

fft_plan fft_create_plan_1d(double complex *in, double complex *out, const unsigned nbits, const int isign, const int FFT_TYPE) {
	fft_plan *plan = malloc(sizeof(fft_plan));
	plan->in = in;
	plan->out = out;
	*(unsigned *) &plan->nbits = nbits;
	*(int *) &plan->isign = isign;
	*(int *) &plan->FFT_TYPE = FFT_TYPE;
	plan->W = malloc(nbits*sizeof(double complex));
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

void fft_danielson_lanczos_mk2(fft_plan plan) {
	memcpy(plan.out, plan.in, (1<<plan.nbits)*sizeof(complex double));
	dfft_bit_reverse(plan.out, plan.nbits);
	int mmax = 1;
	int length = 1<<plan.nbits; 
	int istep;
	double complex ws, temp;
	//printf("\n");
	while (1) {
		if (length <= mmax) break;
		istep = 2*mmax;
		double theta = -M_PI/mmax;
		double complex w = 1.0;
		double complex wp = -2.*sin(0.5*theta)*sin(0.5*theta) + 1.0I*sin(theta);
		//printf("outer loop mmax = %4d istep %4d\n", mmax, istep);
		for (int m = 1; m < mmax+1; m ++) {
			int j = 0;
			ws = w;
			//printf("\tmiddle loop m = %4d\tws = %f, %f\n", m, creal(ws), cimag(ws));
			for (int i = m-1; i < length; i += istep) {
				j = i + mmax;
				temp = ws*plan.out[j];
				plan.out[j] = plan.out[i] - temp;
				plan.out[i] = plan.out[i] + temp;
				//printf("\t\tinner i,j = %d,%d\n", i, j);
			}
			w = w*wp + w;
		} 
		mmax = istep;
	}
}


void fft_danielson_lanczos(fft_plan plan) {
	memcpy(plan.out, plan.in, (1<<plan.nbits)*sizeof(complex double));
	dfft_bit_reverse(plan.out, plan.nbits);
	int n = 1<<plan.nbits;
	for (unsigned s = 0; s < plan.nbits; s++) {
		int m = 1 << (s+1);
		for (int k = 0; k < n; k = k+m) {
			w = 1;
			for (int j = 0; j < m/2; j++) {
				t = w*plan.out[k+j+m/2];
				plan.out[k+j+m/2] = plan.out[k+j] - t;
				plan.out[k+j] = plan.out[k+j] + t;
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
		fft_danielson_lanczos_mk2(plan);
	}
}

void fft_destroy_plan(fft_plan plan) {
	free(plan.tmp);
}

