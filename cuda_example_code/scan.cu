#include <iostream>
#include <cstdio>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdlib>
#include <ctime>
#include <thrust/scan.h>
using namespace std;

  __global__
void prefixScan1(int *in,int *out, int N)
{
  __shared__ int temp[2048];
  int threadId = threadIdx.x;
  int pout = 0, pin = 1;

  //load input into shared memory
  //exclusive scan, so shift right by one and set first element to 0
  temp[threadId] = (threadId > 0) ? in[threadId - 1] : 0;
  __syncthreads();

  for(int offset = 1; offset < N; offset *= 2)
  {
	//swap double buffer
	pout = 1 - pout;
	pin = 1 - pin;

	if(threadId >= offset)
	  temp[pout * N + threadId] = temp[pin * N + threadId] + temp[pin * N + threadId - offset];
	else
	  temp[pout * N + threadId] = temp[pin * N + threadId];
	__syncthreads();
  }

  out[threadId] = temp[pout * N + threadId];
}

  __global__
void prefixScan2(int *in, int *out, int n)
{
  __shared__ int temp[2048];

  int threadId = threadIdx.x;
  int offset = 1;

  //load input into shared memory
  temp[2 * threadId] = in[2 * threadId];
  temp[2 * threadId + 1] = in[2 * threadId + 1];
  __syncthreads();

  for(int d = n/2; d > 0; d /= 2) // build sum in place up the tree
  {
	__syncthreads();
	if(threadId < d)
	{
	  int ai = offset * (2 * threadId + 1) - 1;
	  int bi = offset * (2 * threadId + 2) - 1;
	  temp[bi] += temp[ai];
	}
	offset *= 2;
  }

  if(threadId == 0) // clear the last element
	temp[n-1] = 0;

  for(int d = 1; d < n; d *= 2)
  {
	offset /= 2;
	__syncthreads();

	if(threadId < d)
	{
	  int ai = offset * (2 * threadId + 1) - 1;
	  int bi = offset * (2 * threadId + 2) - 1;
	  int t = temp[ai];
	  temp[ai] = temp[bi];
	  temp[bi] += t;
	}
  }
  __syncthreads();
  out[2 * threadId] = temp[2 * threadId];
  out[2 * threadId + 1] = temp[2 * threadId + 1];
}


int main()
{
  int *h_in,*h_out,*d_in,*d_out;
#ifdef SCAN1
  int N = 1024;
#else
  int N = 2048;
#endif
  
  h_in = (int*)malloc(sizeof(int) * N);
  h_out = (int*)malloc(sizeof(int) * N);

  cudaMalloc(&d_in,sizeof(int) * N);
  cudaMalloc(&d_out,sizeof(int) * N);

  srand(time(NULL));
  for(int i = 0; i < N; ++i)
	h_in[i] = rand() % 5;

  cudaMemcpy(d_in,h_in,sizeof(int) * N, cudaMemcpyHostToDevice);

  dim3 grid(1);
  dim3 block(1024);
  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float ElapsedTime;
  float sumTime = 0;
  for(int i = 0; i < 100; ++i)
  {
	cudaEventRecord(start,0);
#ifdef SCAN1
	prefixScan1<<<grid,block>>>(d_in,d_out,N);
#elif SCAN2
	prefixScan2<<<grid,block>>>(d_in,d_out,N);
#endif
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);
	cudaEventElapsedTime(&ElapsedTime,start,stop);
	sumTime += ElapsedTime;
  }
  printf("average execution time %.3f ms\n ",sumTime / 100);
  cudaMemcpy(h_out,d_out,sizeof(int) * N, cudaMemcpyDeviceToHost);
  thrust::exclusive_scan(h_in,h_in + N,h_in);
  bool success = true;
  for(int i = 0; i < N; ++i)
  {
	if(h_out[i] != h_in[i])
	{
	  success = false;
	  break;
	}
  }

  if(success)
	printf("Success!\n");
  else
	printf("Error!\n");

  free(h_in);
  free(h_out);

  cudaFree(d_in);
  cudaFree(d_out);
  return 0;
}

