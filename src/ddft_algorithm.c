#include "dfft_header.h"
#include <pthread.h>

fft_plan dft_create_plan_1d(double complex *in, double complex *out, const unsigned nbits, const int isign, const unsigned nthreads, const int DFT_TYPE) {
  fft_plan *plan = malloc(sizeof(fft_plan));
  plan->nthreads = nthreads;
  plan->in = in;
  plan->out = out;
  *(unsigned *) &plan->nbits = nbits;
  *(int *) &plan->isign = isign;
  *(int *) &plan->FFT_TYPE = DFT_TYPE;
  plan->W = malloc((1<<nbits)*sizeof(double complex));
  // (serial n^2) algorithm prepare twiddle factors
  double phs = -2.0*M_PI/(1<<nbits);
  for (int s = 0; s < (1<<nbits); s++) {
    plan->W[s] = cexp(-1.0I*phs*s);
  }
  return *plan;
} 

void dft_execute(fft_plan p) {
  long int N = 1<<(p.nbits);
  for (int k = 0; k < N; k++) {
    p.out[k] = p.in[N-1];
    for (int j = 0; j < (1<<p.nbits)-1; j++) {
      p.out[k] = p.out[k]*p.W[k] + p.in[N-2-j];
    }
  }
}

typedef struct {
  fft_plan *p; 
  unsigned thread_id;
} thread_data;

void *thread_work(void* ptr) {
  thread_data *loc = (thread_data *) ptr;
  fft_plan *p = loc->p;
  long int N = 1<<(p->nbits);
  long int Nloc = N/p->nthreads;
  long int istart  = (loc->thread_id)*Nloc;
  long int ifinish = (loc->thread_id+1)*Nloc;
  //printf("Thread %u: works on indx %ld to indx %ld\n", loc->thread_id, istart, ifinish);
  for (int k = istart; k < ifinish; k++) {
    p->out[k] = p->in[N-1];
    for (int j = 0; j < N-1; j++) {
      p->out[k] = p->out[k]*p->W[k] + p->in[N-2-j];
    }
  }
  return NULL;
}

void dft_execute_pthreads(fft_plan p) {
  //long int N = 1<<(p.nbits);
  //long int Nlocal = N/p.nthreads;
  pthread_t *thr_h = malloc(p.nthreads*sizeof(pthread_t));
  thread_data *local_data = malloc(p.nthreads*sizeof(thread_data));
  //printf("\nCreating %d threads\n", p.nthreads);
  for (unsigned i = 0; i < p.nthreads; i++) {
    local_data[i].thread_id = i;
    local_data[i].p = &p;
    pthread_create(&thr_h[i], NULL, thread_work, (void *) &local_data[i]);
  }
  for (unsigned i = 0; i < p.nthreads; i++) {
    pthread_join(thr_h[i], NULL);
  }
}



void dft_destroy_plan(fft_plan plan) {
  free(plan.W);
}
