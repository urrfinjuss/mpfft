#include <stdio.h>
#include <mfft_serial.h>
//#include <math.h>
//#include <complex.h>
#include <time.h>


int main() {
	FILE *fh_log = fopen("mfft_bench.log","a");
	const mpfr_prec_t precision = 128;
	fprintf(fh_log, "# 1. FFT size 2. Avg. time per FFT\n# Precision %u bits\n\n", (unsigned) precision);
	fclose(fh_log);

	mfft_init(precision);
	mfft_version();
	mpfr_printf("Variable Precision with %u bits\n", (unsigned) precision);
	mpfr_printf("Benchmark FFT (Danielson-Lanczos) Performance:\n");
	for (unsigned int nbits = 6; nbits < 11; nbits++) {
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
			//mpfr_div_ui(y, x, 2, MPFR_RNDN);
			//mpfr_tan(y, y, MPFR_RNDN);
			//mpfr_mul(y, y, y, MPFR_RNDN);
			//mpfr_mul_si(y, y, -48, MPFR_RNDN);
			//mpfr_exp(y, y, MPFR_RNDN);
			mpfr_sin(y, x, MPFR_RNDN);
			mpfr_set(f[j].re, y, MPFR_RNDN);
			mpfr_set_ui(f[j].im, 0, MPFR_RNDN);
		}
		// run FFTs
		mfft_plan plan_forward = mfft_create_plan_1d(F, f, nbits, precision, FFT_FORWARD);
		clock_t	begin = clock();
		for (int k = 0; k < 1<<(19-nbits); k++) {
			mfft_execute(plan_forward);
		}
		clock_t	end = clock();
		double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
		//printf("Completed %5d size %4d FFTS in %8.3f secs. Average time per FFT %10.5f ms\n", 1<<(17-nbits), N, time_spent, 1e+3*time_spent/(1<<(17-nbits)));
		printf("Size %5d FFTS in %8.3f secs. Average time per FFT %10.5f ms\n", N, time_spent, 1e+3*time_spent/(1<<(19-nbits)));
		fh_log = fopen("mfft_bench.log","a");
		fprintf(fh_log, "%5d\t%15.10f\n", N, 1e+3*time_spent/(1<<(19-nbits)));
		fclose(fh_log);
		mfft_destroy_plan(&plan_forward);
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
	}

	fh_log = fopen("mfft_bench.log","a");
	fprintf(fh_log, "\n\n");
	fclose(fh_log);


	return 0;
}
