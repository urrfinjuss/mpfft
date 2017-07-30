#include "mfft.h"

static mpc_ptr w;
static mpc_ptr u, t;

void mfft_init(mpfr_prec_t precision) {
	w = init_mpc(precision);
	u = init_mpc(precision);
	t = init_mpc(precision);
	init_bit_operations(precision);
}

mfft_plan mfft_create_plan_1d(mpc_ptr out, mpc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign){
	mfft_plan *p = malloc(sizeof(mfft_plan));
	mpfr_t tmp;

	mpfr_init2(tmp, precision);
	if (out == in) printf("In-Place Transform");
	p->nthreads = 1;
	p->out = out;
	p->in = in;
	p->nbits = nbits;
	p->prec = precision; 
	p->dir = isign;
	mpfr_init2(p->re, precision);
	mpfr_init2(p->im, precision);
	p->W = init_mpc_array(nbits, precision);
	for (unsigned s = 0; s < nbits; s++) {
		mpfr_const_pi(tmp, mode);
		mpfr_div_ui(tmp, tmp, 1 << s, mode);
		mpfr_sin_cos(p->W[s].im, p->W[s].re, tmp, mode);
		mpfr_mul_si (p->W[s].im, p->W[s].im, isign, mode);
	}
	mpfr_clear(tmp);
	return *p;
}

void mfft_destroy_plan(mfft_plan *plan) {
	mpc_clear_array(plan->W, plan->nbits);
	mpfr_clear(plan->re);
	mpfr_clear(plan->im);
	//free(plan);
}


void mfft_execute(mfft_plan plan) {
	mpfr_bit_reverse_copy(plan.out, plan.in, plan.nbits);
	int n = 1<<plan.nbits;

	for (unsigned s = 0; s < plan.nbits; s++) {
		int m = 1 << (s+1);
		//mpfr_printf("Decimation %d: m = %d %.6Re %.6Re\n", s, m, plan.W[s].re,		plan.W[s].im);
		for (int k = 0; k < n; k = k+m) {
			mpfr_set_ui(plan.re, 1, mode);
			mpfr_set_ui(plan.im, 0, mode);
			for (int j = 0; j < m/2; j++) {
				//t = W*plan.out[k+j+m/2]
				mpfr_mul(t->re, plan.im, plan.out[k+j+m/2].im, mode);
				mpfr_fms(t->re, plan.re, plan.out[k+j+m/2].re, t->re, mode);
				mpfr_mul(t->im, plan.im, plan.out[k+j+m/2].re, mode);
				mpfr_fma(t->im, plan.re, plan.out[k+j+m/2].im, t->im, mode);

				//u = plan.out[k+j]
				mpfr_set(u->re, plan.out[k+j].re, mode);
				mpfr_set(u->im, plan.out[k+j].im, mode);

				//A[k+j] = u + t
				mpfr_add(plan.out[k+j].re, u->re, t->re, mode);
				mpfr_add(plan.out[k+j].im, u->im, t->im, mode);

				//A[k+j+m/2] = u - t
				mpfr_sub(plan.out[k+j+m/2].re, u->re, t->re, mode);
				mpfr_sub(plan.out[k+j+m/2].im, u->im, t->im, mode);

				//W = W*Wm;
				mpfr_mul(w->re, plan.im, plan.W[s].im, mode);
				mpfr_fms(w->re, plan.re, plan.W[s].re, w->re, mode);
				mpfr_mul(w->im, plan.im, plan.W[s].re, mode);
				mpfr_fma(w->im, plan.re, plan.W[s].im, w->im, mode);
				mpfr_set(plan.re, w->re, mode);
				mpfr_set(plan.im, w->im, mode);
			}
		}
	}
}


/*
mpc_ptr recursive_fft(mpc_ptr in, long N) {
	if (N == 1) {
		return in;
	} else {
		long m = N/2;
		yt = recursive_fft();
		yb = recursive_fft();
		tw = exp(-2.0I*pi/N).^(0:m-1);
		z = d.*yb;
		y = [yt + z, yt - z];
		return y;
	}
	
}
*/

