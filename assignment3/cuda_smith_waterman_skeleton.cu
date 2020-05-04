#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

#include "cuda_smith_waterman.h"

/*
 *  You can add helper functions and variables as you wish.
 */

__device__
int max_score = 0;

__global__
void kernel(char* d_a /* in */, char* d_b /* in */, 
			int* d_score /* out */, int a_len /* in */, int b_len /* in */, 
			int d /* in; diagonal idx */, int w /* in; width of this diagonal */, int m /* in; mid point of D */
) {
// illegal memory access - allocating too much device memory to CUDA

	// Compute one or more element on a diagonal 'd'
	int i = blockDim.x * blockIdx.x + threadIdx.x;

	// If there are too many threads, only a subset of threads within 'w' will work
	if (i < w) {
		int numCycles = 1;
		// If there are too few threads, thread 'i' will work multiple times cyclically by coalesced access
		int numThreads = gridDim.x * blockDim.x;
		numCycles += (w / numThreads);
		//printf("> (w, d, i) = (%d, %d, %d) has (numThreads, numCycles) = (%d, %d)\n", w, d, i, numThreads, numCycles);
		for (int cycle = 0; cycle < numCycles; cycle++, i += numThreads) {

			// Map from (d, i) -> (x, y) for 'd_a' and 'd_b'
			int delta = max(d - a_len + 1, 0);
			int y = i + delta;
			int x = d - y;

			int current = utils::dev_idx(d, i, b_len);
			/*
			1. i <= m  -> pad for lastlast, lastL; left-parallelogram
			2. i == m+1 -> no pad; forward triangle
			3. i > m+1 -> no pad; right-parallelogram
			*/
			int lastlast, lastL, lastR;
			if (d <= m) {
				lastlast = (i == 0 || i == w-1)	? 0 : d_score[utils::dev_idx(d-2, i-1, b_len)];
				lastL 	 = (i == 0)				? 0 : d_score[utils::dev_idx(d-1, i-1, b_len)];
				lastR 	 = (d == 0) 			? 0 : d_score[utils::dev_idx(d-1, i, b_len)];
			} else {
				lastlast = d_score[utils::dev_idx(d-2, ((d == m+1) ? i : i+1), b_len)];
				lastL 	 = d_score[utils::dev_idx(d-1, i, b_len)];
				lastR	 = d_score[utils::dev_idx(d-1, i+1, b_len)];
			}
			
			d_score[current] = 	max(0,
								max(lastlast + sub_mat(d_a[x], d_b[y]),
								max(lastL - GAP,
									lastR - GAP
								)));
			atomicMax(&max_score, d_score[current]);
			//printf("> (d, i, x, y) = (%d, %d, %d, %d)\n\t(cIdx, L, l, r, d_a[x], d_b[y]) -> (%d, %d, %d, %d, %c, %c)\n\t(score, max_score) = (%d, %d)\n", d, i, x, y, current, lastlast, lastL, lastR, d_a[x], d_b[y], d_score[current], max_score);
		}
	}
}

int smith_waterman(int blocks_per_grid, int threads_per_block, char *a, char *b, int a_len, int b_len) {
	/*
	 *  Please fill in your codes here.
	 */

	// Number of diagonals, Midpoint, Double-max width of a diagonal
	int D = a_len + b_len - 1;
	int m = D / 2;
	int W = 2*m + (D % 2);

	// Allocate device memory
	char *d_a, *d_b;
	// Allocate device-global score 'linear matrix'
	int *d_score;
	cudaMalloc(&d_a, sizeof(char) * a_len);
	cudaMalloc(&d_b, sizeof(char) * b_len);
	cudaMalloc(&d_score, sizeof(int) * (D * b_len));

	// Copy hostToDevice
	cudaMemcpy(d_a, a, sizeof(char) * a_len, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(char) * b_len, cudaMemcpyHostToDevice);
	// Initialise score 'linear matrix' to zeroes
	cudaMemset(d_score, 0, sizeof(int) * (D * b_len));

	// Invoke kernel
	for (int d = 0; d < D; d++) {
		int w = (d <= m) ? (d + 1) : (W - d);
		//printf("(d, w, m) = (%d, %d, %d)\n", d, w, m);
		kernel<<<blocks_per_grid, threads_per_block>>>(d_a, d_b, d_score, a_len, b_len, d, w, m);
	}

	// Return answer
	int answer;
	cudaMemcpyFromSymbol(&answer, max_score, sizeof(int), 0, cudaMemcpyDeviceToHost);
	
	// Free device memory
	cudaFree(d_a);
	cudaFree(d_b);
	cudaFree(d_score);

	return answer;
}
