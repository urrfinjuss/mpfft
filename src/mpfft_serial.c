#include "mpfft_serial.h"

static mpfc_ptr w;
static mpfc_ptr u, t;

void mpfft_init(mpfr_prec_t precision) {
	w = init_mpfc(precision);
	u = init_mpfc(precision);
	t = init_mpfc(precision);
	init_bit_operations(precision);
}

mpfft_plan mpfft_create_plan_1d(mpfc_ptr out, mpfc_ptr in, unsigned nbits, mpfr_prec_t precision, int isign){
	mpfft_plan *p = malloc(sizeof(mpfft_plan));
	mpfr_t tmp;

	mpfr_init2(tmp, precision);
	if (out == in) p->inplace = 1;
	p->out = out;
	p->in = in;
	p->nbits = nbits;
	p->prec = precision; 
	p->dir = isign;
	mpfr_init2(p->re, precision);
	mpfr_init2(p->im, precision);
	p->W = init_mpfc_array(nbits, precision);
	for (unsigned s = 0; s < nbits; s++) {
		mpfr_const_pi(tmp, MODE);
		mpfr_div_ui(tmp, tmp, 1 << s, MODE);
		mpfr_sin_cos(p->W[s].im, p->W[s].re, tmp, MODE);
		mpfr_mul_si (p->W[s].im, p->W[s].im, isign, MODE);
	}
	mpfr_clear(tmp);
	return *p;
}

void mpfft_destroy_plan(mpfft_plan plan) {
	mpfc_clear_array(plan.W, plan.nbits);
	mpfr_clear(plan.re);
	mpfr_clear(plan.im);
}


void mpfft_execute(mpfft_plan plan) {
	if (plan.inplace) mpfr_bit_reverse(plan.out, plan.nbits);
	else mpfr_bit_reverse_copy(plan.out, plan.in, plan.nbits);
	int n = 1<<plan.nbits;

	for (unsigned s = 0; s < plan.nbits; s++) {
		int m = 1 << (s+1);
		//mpfr_printf("Decimation %d: m = %d %.6Re %.6Re\n", s, m, plan.W[s].re,		plan.W[s].im);
		for (int k = 0; k < n; k = k+m) {
			mpfr_set_ui(plan.re, 1, MODE);
			mpfr_set_ui(plan.im, 0, MODE);
			for (int j = 0; j < m/2; j++) {
				//t = W*plan.out[k+j+m/2]
				mpfr_mul(t->re, plan.im, plan.out[k+j+m/2].im, MODE);
				mpfr_fms(t->re, plan.re, plan.out[k+j+m/2].re, t->re, MODE);
				mpfr_mul(t->im, plan.im, plan.out[k+j+m/2].re, MODE);
				mpfr_fma(t->im, plan.re, plan.out[k+j+m/2].im, t->im, MODE);

				//u = plan.out[k+j]
				mpfr_set(u->re, plan.out[k+j].re, MODE);
				mpfr_set(u->im, plan.out[k+j].im, MODE);

				//A[k+j] = u + t
				mpfr_add(plan.out[k+j].re, u->re, t->re, MODE);
				mpfr_add(plan.out[k+j].im, u->im, t->im, MODE);

				//A[k+j+m/2] = u - t
				mpfr_sub(plan.out[k+j+m/2].re, u->re, t->re, MODE);
				mpfr_sub(plan.out[k+j+m/2].im, u->im, t->im, MODE);

				//W = W*Wm;
				mpfr_mul(w->re, plan.im, plan.W[s].im, MODE);
				mpfr_fms(w->re, plan.re, plan.W[s].re, w->re, MODE);
				mpfr_mul(w->im, plan.im, plan.W[s].re, MODE);
				mpfr_fma(w->im, plan.re, plan.W[s].im, w->im, MODE);
				mpfr_set(plan.re, w->re, MODE);
				mpfr_set(plan.im, w->im, MODE);
			}
		}
	}
}


