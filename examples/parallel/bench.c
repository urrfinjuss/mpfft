#include <stdio.h>
#include <mfft.h>
#include <math.h>
#include <time.h>

int main() {
	const mpfr_prec_t precision = 128;
	mfft_pthread_init(precision);
	mfft_version();
	mpfr_printf("Variable Precision with %u bits\n", (unsigned) precision);
	mpfr_printf("Benchmark Threaded FFT (Danielson-Lanczos) Performance:\n");
	unsigned int nbits = 5;
	unsigned int N = 1<<nbits;
	mpc *f = init_mpc_array(N, precision);
	mpc *F = init_mpc_array(N, precision);
	// set initial data
	mpfr_t x, y, scale;
	mpfr_init2(x, precision);
	mpfr_init2(y, precision);
	mpfr_init2(scale, precision);
	mpfr_const_pi(scale, MPFR_RNDN);
	mpfr_div_ui(scale, scale, 1 << (nbits-1), MPFR_RNDN);
	for (int j = 0; j < N; j++) {
		mpfr_mul_si(x, scale, j-N/2, MPFR_RNDN);
		mpfr_sin(y, x, MPFR_RNDN);
		mpfr_set(f[j].re, y, MPFR_RNDN);
		mpfr_set_ui(f[j].im, 0, MPFR_RNDN);
	}
	// run FFTs
	mfft_plan plan_forward = mfft_pthread_create_plan_1d(F, f, 2, nbits, precision, FFT_FORWARD);
	mfft_pthread_execute(plan_forward);
	mfft_pthread_destroy_plan(&plan_forward);
	mpfr_clear(x);
	mpfr_clear(y);
	mpfr_clear(scale);
	char fname[80];
	sprintf(fname, "fft_n%d.txt", 1 << nbits);
	FILE *fhout = fopen(fname, "w");
	for (int j = 0; j < 1<<nbits; j++) {
		mpfr_fprintf(fhout, "%6d\t%36.26Re\t%36.26Re\n", j, F[j].re, F[j].im);
	}
	fclose(fhout);
	mpc_clear_array(f, N);
	mpc_clear_array(F, N);

	return 0;
}
