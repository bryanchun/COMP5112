#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

#include "cuda_smith_waterman.h"

/*
 *  You can add helper functions and variables as you wish.
 */

__global__
void kernel() {

}

int smith_waterman(int blocks_per_grid, int threads_per_block, char *a, char *b, int a_len, int b_len) {
	/*
	 *  Please fill in your codes here.
	 */
	
	// ??Allocate host memory

	// Allocate device memory
	cudaMalloc

	// Copy hostToDevice
	cudaMemcpy(d_data, h_data, sizeof(int)*numElement, 
		   cudaMemcpyHostToDevice);

	// Declare __device__ global score matrix?

	// invoke kernel
	kernel<<<blocks_per_grid, threads_per_block>>>(d_data, numElement);

	// Copy deviceToHost
	cudaMemcpy(h_data, d_data, sizeof(int)*numElement, 
		   cudaMemcpyDeviceToHost);
	
	// Free device memory
	cudaFree(d_data);

	// ??Free host memory

	// Return answer

	// GPUErrChk
	// submat
	// utils::dev_idx

	return -1;
}