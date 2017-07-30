#include "mfft.h"
//#include <omp.h>
#include <pthread.h>
#include <time.h>

static mpc_ptr w;
static mpc_ptr u, t;
static thread_data *data;

void mfft_pthread_init(mpfr_prec_t precision) {
	w = init_mpc(precision);
	u = init_mpc(precision);
	t = init_mpc(precision);
	init_bit_operations(precision);
}

mfft_plan mfft_pthread_create_plan_1d(mpc_ptr out, mpc_ptr in, int nthreads, unsigned nbits, mpfr_prec_t precision, int isign){
	mfft_plan *p = malloc(sizeof(mfft_plan));
	mpfr_t tmp;

	mpfr_init2(tmp, precision);
	if (out == in) printf("In-Place Transform");
	p->nthreads = nthreads;
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

void mfft_pthread_destroy_plan(mfft_plan *plan) {
	mpc_clear_array(plan->W, plan->nbits);
	mpfr_clear(plan->re);
		mpfr_clear(plan->im);
	//free(plan);
}


void mfft_pthread_execute(mfft_plan plan) {
	mpfr_bit_reverse_copy(plan.out, plan.in, plan.nbits);
	int n = 1<<plan.nbits;

	unsigned s = 0;
	int m = 1 << (s+1);
	mpfr_printf("Decimation %d: m = %d %.6Re %.6Re\n", s, m, plan.W[s].re,		plan.W[s].im);
	
	for (int k = 0; k < n; k = k+m) {
		mpfr_set_ui(plan.re, 1, mode);
		mpfr_set_ui(plan.im, 0, mode);
		printf("thread %2d does: %4d and %4d\n", 0, k, k+m/2);
		//t = W*plan.out[k+j+m/2]
		mpfr_mul(t->re, plan.im, plan.out[k+m/2].im, mode);
		mpfr_fms(t->re, plan.re, plan.out[k+m/2].re, t->re, mode);
		mpfr_mul(t->im, plan.im, plan.out[k+m/2].re, mode);
		mpfr_fma(t->im, plan.re, plan.out[k+m/2].im, t->im, mode);
		//u = plan.out[k+j]
		mpfr_set(u->re, plan.out[k].re, mode);
		mpfr_set(u->im, plan.out[k].im, mode);
		//A[k+j] = u + t
		mpfr_add(plan.out[k].re, u->re, t->re, mode);
		mpfr_add(plan.out[k].im, u->im, t->im, mode);
		//A[k+j+m/2] = u - t
		mpfr_sub(plan.out[k+m/2].re, u->re, t->re, mode);
		mpfr_sub(plan.out[k+m/2].im, u->im, t->im, mode);
		//W = W*Wm;
		mpfr_mul(w->re, plan.im, plan.W[s].im, mode);
		mpfr_fms(w->re, plan.re, plan.W[s].re, w->re, mode);
		mpfr_mul(w->im, plan.im, plan.W[s].re, mode);
		mpfr_fma(w->im, plan.re, plan.W[s].im, w->im, mode);
		mpfr_set(plan.re, w->re, mode);
		mpfr_set(plan.im, w->im, mode);
	}

	int nthreads, tid;
	for (unsigned s = 1; s < plan.nbits-1; s++) {
		int m = 1 << (s+1);
		mpfr_printf("Decimation %d: m = %d %.6Re %.6Re\n", s, m, plan.W[s].re,		plan.W[s].im);

		for (int k = 0; k < n; k = k+2*m) {
			mpfr_set_ui(plan.re, 1, mode);
			mpfr_set_ui(plan.im, 0, mode);
			for (int j = 0; j < m/2; j++) {
				printf("thread %2d does: %4d and %4d\n", 0, k+j, k+j+m/2);
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

		for (int k = m; k < n; k = k+2*m) {
			mpfr_set_ui(plan.re, 1, mode);
			mpfr_set_ui(plan.im, 0, mode);
			for (int j = 0; j < m/2; j++) {
				printf("thread %2d does: %4d and %4d\n", 1, k+j, k+j+m/2);
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
	s = plan.nbits-1;
	m = 1 << (s+1);
	mpfr_printf("Decimation %d: m = %d %.6Re %.6Re\n", s, m, plan.W[s].re, plan.W[s].im);
	mpfr_set_ui(plan.re, 1, mode);
	mpfr_set_ui(plan.im, 0, mode);
	for (int j = 0; j < m/2; j++) {
		printf("thread %2d does: %4d and %4d\n", 0, j, j+m/2);
		//t = W*plan.out[k+j+m/2]
		mpfr_mul(t->re, plan.im, plan.out[j+m/2].im, mode);
		mpfr_fms(t->re, plan.re, plan.out[j+m/2].re, t->re, mode);
		mpfr_mul(t->im, plan.im, plan.out[j+m/2].re, mode);
		mpfr_fma(t->im, plan.re, plan.out[j+m/2].im, t->im, mode);

		//u = plan.out[k+j]
		mpfr_set(u->re, plan.out[j].re, mode);
		mpfr_set(u->im, plan.out[j].im, mode);

		//A[k+j] = u + t
		mpfr_add(plan.out[j].re, u->re, t->re, mode);
		mpfr_add(plan.out[j].im, u->im, t->im, mode);

		//A[k+j+m/2] = u - t
		mpfr_sub(plan.out[j+m/2].re, u->re, t->re, mode);
		mpfr_sub(plan.out[j+m/2].im, u->im, t->im, mode);

		//W = W*Wm;
		mpfr_mul(w->re, plan.im, plan.W[s].im, mode);
		mpfr_fms(w->re, plan.re, plan.W[s].re, w->re, mode);
		mpfr_mul(w->im, plan.im, plan.W[s].re, mode);
		mpfr_fma(w->im, plan.re, plan.W[s].im, w->im, mode);
		mpfr_set(plan.re, w->re, mode);
		mpfr_set(plan.im, w->im, mode);
	}
	
}

void *thread(void *arg) {
	thread_data *local = (thread_data *) arg;
	mpfr_prec_t precision = local->precision;

	mpfr_init2(local->sum, precision);
	mpfr_set_si(local->sum, 0, mode);

	int tid = local->tid;
	int stride = local->stride;
	printf("Thread %d: Start Index %4d\tEnd Index %4d\n", tid, tid*stride, (tid+1)*stride);
	for (int j = tid*stride; j < (tid+1)*stride; j++) {
		//local->sum += j*dx*(1.0 - j*dx);
		mpfr_add(local->sum, local->sum, local->f[j].im, mode);
		//local->sum += sin(2.*M_PI*j*dx);
	}
	pthread_exit(NULL);
}

void mfft_pthread_example(const int nthreads, mpfr_prec_t precision) {
	pthread_t     tid[nthreads];
	data = malloc(nthreads*sizeof(thread_data));
	int nbits = 20;
	int N = 1<<nbits;
	int stride = N/nthreads;

	printf("N = %d\nnthreads = %d\nstride = %d\n\n", N, nthreads, stride);
	// Set Global Data
	mpfr_t x, dx;
	mpc_ptr F = init_mpc_array(N, precision);
	mpfr_init2(x, precision);
	mpfr_init2(dx, precision);
	mpfr_const_pi(dx, mode);
	mpfr_div_si(dx, dx, N, mode);
	mpfr_mul_ui(dx, dx, 2, mode);
	for (int j = 0; j < N; j++) {
		mpfr_mul_si(x, dx, j, mode);
		mpfr_sin_cos(F[j].im, F[j].re, x, mode);
	}  
	// End Set Global
	// Set Local Data & Run Multithreaded
	int i;
	for (i = 0; i < nthreads; i++) {
		data[i].tid = i;
		data[i].stride = stride;
		data[i].precision = precision;
		mpfr_init2(data[i].dx, precision);
		mpfr_set(data[i].dx, dx, mode);
		data[i].f = F;
	}
	struct timespec start, finish;
	double elapsed;
	clock_gettime(CLOCK_MONOTONIC, &start);

	for (i = 0; i < nthreads; i++) {
		if (pthread_create(&tid[i], NULL, thread, (void *) &data[i])) {
			printf("pthread_create() failed\n");
			exit(-1);
		}; 
	}

	mpfr_t full;
	mpfr_init2(full, precision);
	mpfr_set_si(full, 0, mode);
	for (i = 0; i < nthreads; i++) {
		pthread_join(tid[i], NULL); 
		mpfr_fma(full, data[i].sum, dx, full, mode);
	}
	clock_gettime(CLOCK_MONOTONIC, &finish);
	elapsed = (finish.tv_sec - start.tv_sec);
	elapsed += (finish.tv_nsec - start.tv_nsec) / 1000000000.0;
	//mpfr_printf("Full Sum is %.32Re\n\n", full);
	printf("Size %5d MPFR Array %10.5f ms on %d threads\n", N, 1e+3*elapsed, nthreads);
}
