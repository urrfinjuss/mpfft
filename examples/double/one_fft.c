#include <stdio.h>
#include <base_fft.h>
#include <math.h>
#include <complex.h>
#include <time.h>


int main() {
	mfft_version();
	int nbits = 4;
	int N = 1<<nbits;
	double complex *f = malloc(N*sizeof(double complex));
	double complex *F = malloc(N*sizeof(double complex));
	// set initial data
	double scale = 2*M_PI/N;
	double x;
	for (int j = 0; j < N; j++) {
		x = (j-N/2)*scale;
		f[j] = sin(x) + 0.0I;
		printf("%e,\t", creal(f[j]));
	}
	printf("\n");
	// run FFTs
	fft_plan plan_forward = fft_create_plan_1d(f, F, N, FFT_FORWARD);
	fft_execute(plan_forward);
	fft_destroy_plan(plan_forward);

	char fname[80];
	sprintf(fname, "fft_n%d.txt", 1 << nbits);
	FILE *fhout = fopen(fname, "w");
	for (int j = 0; j < 1<<nbits; j++) {
		fprintf(fhout, "%6d\t%22.15e\t%22.15e\n", j, creal(F[j]), cimag(F[j]));
	}
	fclose(fhout);
	free(f);
	free(F);

	return 0;
}
