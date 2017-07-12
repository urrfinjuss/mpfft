#include "mfft.h"

naive_plan init_mfft_naive_plan(mpc_ptr out, mpc_ptr in, unsigned N, mpfr_prec_t precision, int isign) {
	naive_plan	*plan = malloc(sizeof(naive_plan));
	if (out == in) 	printf("In-Place Transform %u points\n", N);
	else 						printf("Out-of-Place Transform %u points\n", N);

	mpfr_init2(plan->tmp1, precision);
	mpfr_init2(plan->tmp2, precision);
	plan->N = N;
	plan->w = init_mpc_array(N, precision);
	plan->in = in;
	plan->out = out;
	mpfr_set_ui(plan->w[0].re, 1, mode);
	mpfr_set_ui(plan->w[0].im, 0, mode);

	mpfr_const_pi(plan->tmp1, precision);
	mpfr_mul_ui  (plan->tmp1, plan->tmp1, 2, mode);
	mpfr_div_ui  (plan->tmp1, plan->tmp1, N, mode);

	mpfr_sin_cos(plan->w[1].im, plan->w[1].re, plan->tmp1, mode);
	mpfr_mul_si (plan->w[1].im, plan->w[1].im, isign, mode);

	for (unsigned k = 1; k < N-1; k++) {
		mpfr_mul(plan->w[k+1].re, plan->w[k].im, plan->w[1].im, mode);
		mpfr_fms(plan->w[k+1].re, plan->w[k].re, plan->w[1].re, plan->w[k+1].re, mode);
		mpfr_mul(plan->w[k+1].im, plan->w[k].re, plan->w[1].im, mode);
		mpfr_fma(plan->w[k+1].im, plan->w[k].im, plan->w[1].re, plan->w[k+1].im, mode);
	}
	return *plan;
}

void mfft_naive_execute(naive_plan p) {
	unsigned n = p.N;
	for (unsigned k = 0; k < n; k++) {
		mpfr_mul(p.tmp1, p.in[n-1].im, p.w[k].im, mode);
		mpfr_fms(p.tmp1, p.in[n-1].re, p.w[k].re, p.tmp1, mode); 

		mpfr_mul(p.tmp2, p.in[n-1].re, p.w[k].im, mode);
		mpfr_fma(p.tmp2, p.in[n-1].im, p.w[k].re, p.tmp2, mode); 

		mpfr_set(p.out[k].re, p.tmp1, mode);
		mpfr_set(p.out[k].im, p.tmp2, mode);
		for (unsigned j = 1; j < n-1; j++) {
			mpfr_add(p.out[k].re, p.out[k].re, p.in[n-j-1].re, mode);
			mpfr_add(p.out[k].im, p.out[k].im, p.in[n-j-1].im, mode);

			mpfr_mul(p.tmp1, p.out[k].im, p.w[k].im, mode);
			mpfr_fms(p.tmp1, p.out[k].re, p.w[k].re, p.tmp1, mode); 

			mpfr_mul(p.tmp2, p.out[k].re, p.w[k].im, mode);
			mpfr_fma(p.tmp2, p.out[k].im, p.w[k].re, p.tmp2, mode); 
			mpfr_set(p.out[k].re, p.tmp1, mode);
			mpfr_set(p.out[k].im, p.tmp2, mode);
		}
		mpfr_add(p.out[k].re, p.out[k].re, p.in[0].re, mode);
		mpfr_add(p.out[k].im, p.out[k].im, p.in[0].im, mode);
	}
}
