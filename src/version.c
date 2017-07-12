#include "mfft.h"

void mfft_version() {
	printf("MPFR FFT version %u.%u\n", MFFT_MAJOR_VERSION, MFFT_MINOR_VERSION);
	printf("MPFR Version %s\n", mpfr_get_version());
}



