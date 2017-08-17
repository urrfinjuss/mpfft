
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math_constants.h"

#include"cuComplex.h"
#include <stdlib.h>
#include <stdio.h>


//A wrapper function to call cuda safely (A ton of error checks)
cudaError_t cudaDDFT(unsigned nbits, unsigned nthreads, cuDoubleComplex * out, cuDoubleComplex * in);

__global__ void ddftKernel(unsigned nbits, unsigned nthreads, cuDoubleComplex* out, cuDoubleComplex* in, cuDoubleComplex* W)
{
	long int N = 1 << nbits;
	long int Nloc = N / nthreads;
	long int tid = threadIdx.x;
	long int istart = tid*Nloc;
	long int ifinish = (tid+1)*Nloc;
	for (int k = istart; k < ifinish; k++) {
		out[k] = in[N - 1];
		for (int j = 0; j < N - 1; j++) {
			out[k].x = out[k].x * W[k].x - out[k].y * W[k].y + in[N - 2 - j].x;
			out[k].y = out[k].y * W[k].x + out[k].x * W[k].y + in[N - 2 - j].y;
		}
	}
}

int main()
{
	unsigned nthreads = 16;
	unsigned nbits = 10;

	//Create test
	cuDoubleComplex* in = (cuDoubleComplex*)malloc(1 << nbits * sizeof(cuDoubleComplex));
	cuDoubleComplex* out = (cuDoubleComplex*)malloc(1 << nbits * sizeof(cuDoubleComplex));

	FILE *fh = fopen("init.txt", "w");
	for (unsigned i = 0; i < 1 << nbits; ++i)
	{
		in[i].x = sin(2. * i * CUDART_PI / (1 << nbits));
		in[i].y = 0;
		fprintf(fh, "%19.12e\t%19.12e\t%19.12e\n", 2.*CUDART_PI*i / (1 << nbits), in[i].x, in[i].y);
	}
	fclose(fh);


	//Measure time


    cudaError_t cudaStatus = cudaDDFT(nbits, nthreads, out, in);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addWithCuda failed!");
        return 1;
    }

	fh = fopen("output.txt", "w");
	for (unsigned i = 0; i < 1 << nbits; ++i)
	{
		fprintf(fh, "%19.12e\t%19.12e\t%19.12e\n", 2.*CUDART_PI*i / (1 << nbits), out[i].x, out[i].y);
	}
	fclose(fh);

    // cudaDeviceReset must be called before exiting in order for profiling and
    // tracing tools such as Nsight and Visual Profiler to show complete traces.
    cudaStatus = cudaDeviceReset();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceReset failed!");
        return 1;
    }

    return 0;
}

// Helper function for using CUDA to compute discrete Fourier transform in parallel.
cudaError_t cudaDDFT(unsigned nbits, unsigned nthreads, cuDoubleComplex * out, cuDoubleComplex * in)
{
	cuDoubleComplex* dev_in;
	cuDoubleComplex* dev_out;
	cuDoubleComplex* dev_W;
	cuDoubleComplex* W;

	W = (cuDoubleComplex*)malloc((1 << nbits) * sizeof(cuDoubleComplex));
	// (serial n^2) algorithm prepare twiddle factors
	double phs = -2.0*CUDART_PI / (1 << nbits);
	for (int s = 0; s < (1 << nbits); s++) 
	{
		W[s].x = cos(-phs*s);
		W[s].y = sin(-phs*s);
	}
	cudaError_t cudaStatus;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaStatus = cudaSetDevice(0);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaSetDevice failed!  Do you have a CUDA-capable GPU installed?");
        goto Error;
    }

    // Allocate GPU buffers for two vectors (one input, one output)    .
    cudaStatus = cudaMalloc((void**)&dev_out, (1 << nbits) * sizeof(cuDoubleComplex));
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMalloc failed!");
        goto Error;
    }

	cudaStatus = cudaMalloc((void**)&dev_in, (1 << nbits) * sizeof(cuDoubleComplex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}

	cudaStatus = cudaMalloc((void**)&dev_W, (1 << nbits) * sizeof(cuDoubleComplex));
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMalloc failed!");
		goto Error;
	}


    // Copy input vectors from host memory to GPU buffers.
    cudaStatus = cudaMemcpy(dev_in, in, (1<<nbits) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

	cudaStatus = cudaMemcpy(dev_W, W, (1 << nbits) * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice);
	if (cudaStatus != cudaSuccess) {
		fprintf(stderr, "cudaMemcpy failed!");
		goto Error;
	}

    // Launch a kernel on the GPU with one thread for each element.
    ddftKernel<<<1, nthreads>>>(nbits, nthreads, dev_out, dev_in, dev_W);

    // Check for any errors launching the kernel
    cudaStatus = cudaGetLastError();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "addKernel launch failed: %s\n", cudaGetErrorString(cudaStatus));
        goto Error;
    }
    
    // cudaDeviceSynchronize waits for the kernel to finish, and returns
    // any errors encountered during the launch.
    cudaStatus = cudaDeviceSynchronize();
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaDeviceSynchronize returned error code %d after launching addKernel!\n", cudaStatus);
        goto Error;
    }

    // Copy output vector from GPU buffer to host memory.
    cudaStatus = cudaMemcpy(out, dev_out, (1 << nbits) * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost);
    if (cudaStatus != cudaSuccess) {
        fprintf(stderr, "cudaMemcpy failed!");
        goto Error;
    }

Error:
    cudaFree(dev_out);
    cudaFree(dev_in);
    cudaFree(dev_W);
    
    return cudaStatus;
}
