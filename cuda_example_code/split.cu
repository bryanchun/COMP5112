/*
   Parallel split kernels
 */
#include <iostream>
#include <cstdio>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdlib>
#include <ctime>
#include <helper_cuda.h>
#include <thrust/scan.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
using namespace std;

__device__
int getPartID(int element)
{
  return element % 2;
}

__global__ 
void mapPart(int *d_R,int *d_pidArray,int r_len)
{
  int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  int threadNumber = blockDim.x * gridDim.x;
  while(threadId < r_len)
  {
	d_pidArray[threadId] = getPartID(d_R[threadId]);
	threadId += threadNumber;
  }
}

__global__
void count_Hist(int *d_Hist,int *d_pidArray,int r_len, int numPart)
{
  __shared__ int s_Hist[2048];
  int threadId = blockIdx.x * blockDim.x + threadIdx.x;
  int threadNumber = blockDim.x * gridDim.x;
  int offset = threadIdx.x * numPart;
  
  for(int i = 0; i < numPart; ++i)
	s_Hist[i + offset] = 0;
  
  for(int i = threadId; i < r_len; i += threadNumber)
	s_Hist[offset + d_pidArray[i]]++;

  for(int i = 0; i < numPart; ++i)
	d_Hist[i * threadNumber + threadId] = s_Hist[offset + i];
}

__global__
void write_Hist(int d_pidArray[],int d_psSum[],int d_loc[],int numPart,int r_len)
{
  __shared__ int s_psSum[2048];
  int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  int threadNumber = gridDim.x * blockDim.x;
  int offset = threadIdx.x * numPart;

  for(int i = 0; i < numPart; ++i)
	s_psSum[i + offset] = d_psSum[threadId + i * threadNumber];

  for(int i = threadId; i < r_len; i += threadNumber)
  {
	int pid = d_pidArray[i];
	d_loc[i] = s_psSum[pid + offset];
	s_psSum[pid+ offset]++;
  }
}

__global__
void scatter(int d_Rin[],int d_Rout[],int d_loc[],int r_len)
{
  int threadId = threadIdx.x + blockIdx.x * blockDim.x;
  int threadNumber = blockDim.x * gridDim.x;

  while(threadId < r_len)
  {
	d_Rout[d_loc[threadId]] = d_Rin[threadId];
	threadId += threadNumber;
  }
}

int main()
{
  int *d_R,*d_result,*d_pidArray,*d_loc,*d_Hist,*d_psSum;
  int *h_R,*h_result;
  int r_len = 1024 * 1024;
  int numPart = 2; //even or odd
  dim3 grid(256);
  dim3 block(512);
  int numThread = grid.x * block.x;
  int Hist_len = numThread * numPart;

  h_R = (int*) malloc(sizeof(int) * r_len);
  h_result = (int*) malloc(sizeof(int) * r_len);
  checkCudaErrors(cudaMalloc(&d_R,sizeof(int) * r_len));
  checkCudaErrors(cudaMalloc(&d_result,sizeof(int) * r_len));
  checkCudaErrors(cudaMalloc(&d_pidArray,sizeof(int) * r_len));
  checkCudaErrors(cudaMalloc(&d_loc,sizeof(int) * r_len));
  checkCudaErrors(cudaMalloc(&d_Hist,sizeof(int) * Hist_len));
  checkCudaErrors(cudaMalloc(&d_psSum,sizeof(int) * Hist_len));

  
  srand(0);
  for(int i = 0; i < r_len; ++i)
	h_R[i] = rand() % 2;
  cudaMemcpy(d_R,h_R,sizeof(int) * r_len, cudaMemcpyHostToDevice);

  mapPart<<<grid,block>>>(d_R,d_pidArray,r_len);
  count_Hist<<<grid,block>>>(d_Hist,d_pidArray,r_len,numPart);

  thrust::device_ptr<int> dev_Hist(d_Hist);
  thrust::device_ptr<int> dev_psSum(d_psSum);
  thrust::exclusive_scan(dev_Hist,dev_Hist + Hist_len,dev_psSum);

  write_Hist<<<grid,block>>>(d_pidArray,d_psSum,d_loc,numPart,r_len);
  scatter<<<grid,block>>>(d_R,d_result,d_loc,r_len);
  
  cudaMemcpy(h_result,d_result,sizeof(int) * r_len,cudaMemcpyDeviceToHost);
  
  freopen("before.txt","w",stdout);
  for(int i = 0; i < r_len; ++i)
	printf("%d %d\n",i,h_R[i]);
  fclose(stdout);
  freopen("/dev/tty/","w",stdout);

  freopen("after.txt","w",stdout);
  for(int i = 0; i < r_len; ++i)
	printf("%d %d\n",i,h_result[i]);
  fclose(stdout);
  freopen("/dev/tty","w",stdout);

  free(h_R);
  free(h_result);
  checkCudaErrors(cudaFree(d_R));
  checkCudaErrors(cudaFree(d_result));
  checkCudaErrors(cudaFree(d_pidArray));
  checkCudaErrors(cudaFree(d_loc));
  checkCudaErrors(cudaFree(d_Hist));
  checkCudaErrors(cudaFree(d_psSum));
 
  return 0;
}
