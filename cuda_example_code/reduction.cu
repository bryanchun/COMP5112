/*
   Parallel reduction kernels
 */
#include <iostream>
#include <cstdio>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdlib>
using namespace std;

/* In this version only those threads whose IDs are multiples of the stride work 
*/
  __global__ void
reduce0(int *g_idata, int *g_odata, unsigned int n)
{
  __shared__ int sdata[1024];

  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

  sdata[tid] = (i < n) ? g_idata[i] : 0;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=1; s < blockDim.x; s *= 2)
  {
        //check if current thread should work in this iteration
        if ((tid % (2*s)) == 0)
          sdata[tid] += sdata[tid + s];
        __syncthreads();
  }
  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}


/* 
   In this version working threads are consecutive 
 */
  __global__ 
void reduce1(int *g_idata, int *g_odata, unsigned int n)
{
  __shared__ int sdata[1024];
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;

  sdata[tid] = (i < n) ? g_idata[i] : 0;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=1; s < blockDim.x; s *= 2)
  {
        int index = 2 * s * tid;
        if (index < blockDim.x)
          sdata[index] += sdata[index + s];
        __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

/*
   In this version consecutive threads perform coalesced access
 */
  __global__ void
reduce2(int *g_idata, int *g_odata, unsigned int n)
{

  __shared__ int sdata[1024];
  // load shared mem
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*blockDim.x + threadIdx.x;
  sdata[tid] = (i < n) ? g_idata[i] : 0;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
        if (tid < s)
        {
          sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = sdata[0];
}

/*
   this version uses n/2 threads --
   it performs the first level of reduction when reading from global memory.
 */
  __global__ void
reduce3(int *g_idata, int *g_odata, unsigned int n)
{
  __shared__ int sdata[1024];

  // perform first level of reduction,
  // reading from global memory, writing to shared memory
  unsigned int tid = threadIdx.x;
  unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;

  int mySum = (i < n) ? g_idata[i] : 0;

  if (i + blockDim.x < n)
        mySum += g_idata[i+blockDim.x];

  sdata[tid] = mySum;
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=blockDim.x/2; s>0; s>>=1)
  {
        if (tid < s)
        {
          sdata[tid] = mySum = mySum + sdata[tid + s];
        }

        __syncthreads();
  }

  // write result for this block to global mem
  if (tid == 0) g_odata[blockIdx.x] = mySum;
}

int main()
{
  int N = 64 * 1024 * 1024;
  int dataSize = sizeof(int) * N;
#ifdef REDUCE3
  dim3 grid(32 * 1024);
#else
  dim3 grid(64 * 1024);
#endif
  dim3 block(1024);
  int *in,*out;
  int *d_in,*d_out;
  //grid.x is the number of blocks
  in = (int*)malloc(sizeof(int) * N);
  out = (int*)malloc(sizeof(int) * grid.x);

  cudaMalloc(&d_in,sizeof(int) * N);
  cudaMalloc(&d_out,sizeof(int) * grid.x);

  for(int i = 0; i < N; ++i)
        in[i] = 1;

  cudaMemcpy(d_in,in,sizeof(int) * N, cudaMemcpyHostToDevice);

  cudaEvent_t start,stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  float ElapsedTime;
  float sumTime = 0.0;

  for(int i = 0; i < 100; ++i)
  {
        cudaEventRecord(start,0);
#ifdef REDUCE0
        reduce0<<<grid,block>>>(d_in,d_out,N);
#elif REDUCE1
        reduce1<<<grid,block>>>(d_in,d_out,N);
#elif REDUCE2
        reduce2<<<grid,block>>>(d_in,d_out,N);
#else
        reduce3<<<grid,block>>>(d_in,d_out,N);
#endif
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&ElapsedTime,start,stop);
        sumTime += ElapsedTime;
  }

  double throughPut = dataSize / (sumTime / 100 * 1e-3);
  printf("average execution time %.3f ms\nThroughput %.4lf GB/s\n ",sumTime / 100,throughPut * 1e-9);
  
  cudaMemcpy(out,d_out,sizeof(int) * grid.x,cudaMemcpyDeviceToHost);
  int sum = 0;
  for(int i = 0; i < grid.x; ++i)
        sum += out[i];
  if(sum == N)
        printf("Success! sum is %d, out[0] is %d\n", sum, out[0]);
  else
        printf("Error! sum is %d, out[0] is %d\n", sum, out[0]);

  free(in);
  free(out);

  cudaFree(d_in);
  cudaFree(d_out);
  return 0;
}
