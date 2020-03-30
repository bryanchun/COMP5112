/*
 Parallel Radix Sort
 */
#include <iostream>
#include <algorithm>
#include <cstdio>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <cstdlib>
#include <ctime>
#include <sys/time.h>
#include <helper_cuda.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>
using namespace std;

__device__
int getPartID(int element, int mask, int shiftBits) {
	element >>= shiftBits;
	return element & mask;
}

__global__
void mapPart(int *d_R, int *d_pidArray, int r_len, int shiftBits, int mask) {
	int threadId = blockIdx.x * blockDim.x + threadIdx.x;
	int threadNumber = blockDim.x * gridDim.x;
	while (threadId < r_len) {
		d_pidArray[threadId] = getPartID(d_R[threadId], mask, shiftBits);
		threadId += threadNumber;
	}
}

__global__
void count_Hist(int *d_Hist, int *d_pidArray, int r_len, int numPart) {
	__shared__ int s_Hist[512 * 16];
	int threadId = blockIdx.x * blockDim.x + threadIdx.x;
	int threadNumber = blockDim.x * gridDim.x;
	int elePerThread = r_len / threadNumber;
	int offset = threadIdx.x * numPart;

	for (int i = 0; i < numPart; ++i)
		s_Hist[i + offset] = 0;

	for (int i = 0; i < elePerThread; ++i)
		s_Hist[offset + d_pidArray[i + threadId * elePerThread]]++;

	for (int i = 0; i < numPart; ++i)
		d_Hist[i * threadNumber + threadId] = s_Hist[offset + i];
}

__global__
void write_Hist(int d_pidArray[], int d_psSum[], int d_loc[], int numPart,
		int r_len) {
	__shared__ int s_psSum[512 * 16];
	int threadId = threadIdx.x + blockIdx.x * blockDim.x;
	int threadNumber = gridDim.x * blockDim.x;
	int elePerThread = r_len / threadNumber;
	int offset = threadIdx.x * numPart;

	for (int i = 0; i < numPart; ++i)
		s_psSum[i + offset] = d_psSum[threadId + i * threadNumber];

	for (int i = 0; i < elePerThread; ++i) {
		int pid = d_pidArray[i + threadId * elePerThread];
		d_loc[i + threadId * elePerThread] = s_psSum[pid + offset];
		s_psSum[pid + offset]++;
	}
}

__global__
void scatter(int d_Rin[], int d_Rout[], int d_loc[], int r_len) {
	int threadId = threadIdx.x + blockIdx.x * blockDim.x;
	int threadNumber = blockDim.x * gridDim.x;

	while (threadId < r_len) {
		d_Rout[d_loc[threadId]] = d_Rin[threadId];
		threadId += threadNumber;
	}
}

double diffTime(timeval start, timeval end) {
	return (end.tv_sec - start.tv_sec) * 1000
			+ (end.tv_usec - start.tv_usec) * 0.001;
}
int main() {
	struct timeval cpu_start, cpu_end;
	struct timeval gpu_start, gpu_end;
	int *d_R, *d_result, *d_pidArray, *d_loc, *d_Hist, *d_psSum, *d_resultTmp;
	int *h_R, *h_result;
	int r_len = 1024 * 1024;

	int numBits = 16;
	int bitsPerPass = 4;
	int numPass = numBits / bitsPerPass;
	int mask = (1 << bitsPerPass) - 1;

	int numPart = 1 << bitsPerPass;
	dim3 grid(256);
	dim3 block(512);
	int numThread = grid.x * block.x;
	int Hist_len = numThread * numPart;

	h_R = (int*) malloc(sizeof(int) * r_len);
	h_result = (int*) malloc(sizeof(int) * r_len);
	checkCudaErrors(cudaMalloc(&d_R, sizeof(int) * r_len));
	checkCudaErrors(cudaMalloc(&d_result, sizeof(int) * r_len));
	checkCudaErrors(cudaMalloc(&d_resultTmp, sizeof(int) * r_len));
	checkCudaErrors(cudaMalloc(&d_pidArray, sizeof(int) * r_len));
	checkCudaErrors(cudaMalloc(&d_loc, sizeof(int) * r_len));
	checkCudaErrors(cudaMalloc(&d_Hist, sizeof(int) * Hist_len));
	checkCudaErrors(cudaMalloc(&d_psSum, sizeof(int) * Hist_len));

	thrust::device_ptr<int> dev_Hist(d_Hist);
	thrust::device_ptr<int> dev_psSum(d_psSum);

	srand(time(NULL));
	for (int i = 0; i < r_len; ++i)
		h_R[i] = rand() % (1 << numBits);

	gettimeofday(&gpu_start, NULL);
	cudaDeviceSynchronize();
	cudaMemcpy(d_R, h_R, sizeof(int) * r_len, cudaMemcpyHostToDevice);
	for (int passID = 0; passID < numPass; ++passID) {
		int shiftBits = passID * bitsPerPass;
		mapPart<<<grid,block>>>(d_R,d_pidArray,r_len,shiftBits,mask);
		count_Hist<<<grid,block>>>(d_Hist,d_pidArray,r_len,numPart);
		thrust::exclusive_scan(dev_Hist, dev_Hist + Hist_len, dev_psSum);
		write_Hist<<<grid,block>>>(d_pidArray,d_psSum,d_loc,numPart,r_len);
		scatter<<<grid,block>>>(d_R,d_result,d_loc,r_len);
		cudaMemcpy(d_R, d_result, sizeof(int) * r_len,
				cudaMemcpyDeviceToDevice);
	}
	cudaMemcpy(h_result, d_result, sizeof(int) * r_len, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	gettimeofday(&gpu_end, NULL);
	printf("GPU execution time %.3fms\n", diffTime(gpu_start, gpu_end));

	gettimeofday(&cpu_start, NULL);
	sort(h_R, h_R + r_len);
	gettimeofday(&cpu_end, NULL);
	printf("CPU execution time %.3fms\n", diffTime(cpu_start, cpu_end));

	bool success = true;
	for (int i = 0; i < r_len; ++i) {
		if (h_R[i] != h_result[i]) {
			success = false;
			break;
		}
	}
	if (success)
		printf("Success!\n");
	else
		printf("Error!\n");

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

