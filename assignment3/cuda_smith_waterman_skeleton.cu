#include <iostream>
#include <string>
#include <algorithm>
#include <chrono>

using namespace std;

#include "cuda_smith_waterman.h"

/*
 *  You can add helper functions and variables as you wish.
 */

__const__
int d_tileHeight = 10000;
const int h_tileHeight = 10000;
const int MAX_SEQ_SIZE = 20000;

__device__
int max_score = 0;

// +2 for offsetting tile's top 2 rows
inline __device__
int lin_idx(int x, int y, int n) {
	return utils::dev_idx(x+2, y, n);
}

__global__
void kernel(
	char* d_a /* in */, char* d_b /* in */, int a_len /* in */, int b_len /* in */, 
	int* d_score /* out */,
	int d /* in; diagonal idx */, int w /* in; width of this diagonal */, int D /* in; number of diagonals */
) {
// Incorrect input3.txt?

	// Compute one or more element on a diagonal 'd'
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	int d_tiled = d % d_tileHeight;

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

			int current = lin_idx(d_tiled, i, b_len);
			/*
			1. d < a_len  -> only pad for lastlast, lastL outside of leftmost bound; left-parallelogram
			2. d == a_len -> no pad; forward triangle
			3. d > a_len  -> no pad; right-parallelogram
			*/
			int lastlast, lastL, lastR;
			if (d < a_len) {
				lastlast = (i == 0) ? 0 : d_score[lin_idx(d_tiled-2, i-1, b_len)];
				lastL 	 = (i == 0) ? 0 : d_score[lin_idx(d_tiled-1, i-1, b_len)];
				lastR 	 = d_score[lin_idx(d_tiled-1, i, b_len)];
			} else {
				lastlast = d_score[lin_idx(d_tiled-2, ((d == a_len) ? i : i+1), b_len)];
				lastL 	 = d_score[lin_idx(d_tiled-1, i, b_len)];
				lastR	 = d_score[lin_idx(d_tiled-1, i+1, b_len)];
				//printf("* (d, i, x, y) = (%d, %d, %d, %d): L of (%d, %d), l of (%d, %d), r of (%d, %d)\n", d, i, x, y, d-2, ((d == a_len) ? i : i+1), d-1, i, d-1, i+1);
			}
			
			d_score[current] =  max(0,
						max(lastlast + sub_mat(d_a[x], d_b[y]),
						max(lastL - GAP,
							lastR - GAP
						)));
			atomicMax(&max_score, d_score[current]);
			//printf("> (d, d_tiled, i, x, y) = (%d, %d, %d, %d, %d)\n\t(cIdx, L, l, r, d_a[x], d_b[y]) -> (%d, %d, %d, %d, %c, %c)\n\t(score, max_score) = (%d, %d)\n", d, d_tiled, i, x, y, current, lastlast, lastL, lastR, d_a[x], d_b[y], d_score[current], max_score);
		}
	}
}

__global__
void nextTile(int* d_score, int d_score_height, int b_len) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	// If there are too many threads, only a subset of threads within 'b_len' will work
	if (i < b_len) {
		int numCycles = 1;
		// If there are too few threads, thread 'i' will work multiple times cyclically by coalesced access
		int numThreads = gridDim.x * blockDim.x;
		numCycles += (b_len / numThreads);
		for (int cycle = 0; cycle < numCycles; cycle++, i += numThreads) {
			d_score[utils::dev_idx(0, i, b_len)] = d_score[utils::dev_idx(d_score_height-2, i, b_len)];
			d_score[utils::dev_idx(1, i, b_len)] = d_score[utils::dev_idx(d_score_height-1, i, b_len)];
			for (int y = 2; y < d_score_height; y++) {
				d_score[utils::dev_idx(y, i, b_len)] = 0;
			}
		}
	}
}

__global__
void debug_xy(char* d_a, char* d_b, int* d_score, int a_len, int b_len) {
	printf("\t");
	for (int i = 0; i < b_len; i++) {
        printf("%c\t", d_b[i]);
	}
	printf("\n");

	for (int i = 1; i <= a_len; i++) {
		printf("%c\t", d_a[i-1]);
		for (int j = 1; j <= b_len; j++) {
			int d = (i-1) + (j-1);
			int i_prime = (j-1) - max(d - a_len + 1, 0);
			printf("%d\t", d_score[lin_idx(d, i_prime, b_len)]);
		}
		printf("\n");
	}
	printf("\n");
}

// __global__
// void debug_di(int* d_score, int d_score_height, int b_len) {
// 	for (int i = 0; i < d_score_height; i++) {
// 		for (int j = 0; j < b_len; j++) {
// 			printf("%d\t", d_score[utils::dev_idx(j, i, b_len)]);
// 		}
// 		printf("\n");
// 	}
// 	printf("\n");
// }

int smith_waterman(int blocks_per_grid, int threads_per_block, char *a, char *b, int a_len, int b_len) {
	/*
	 *  Please fill in your codes here.
	 */
	
	// Number of diagonals, Double-max width of a diagonal
	int D = a_len + b_len - 1;
	int W = 2*a_len + (D % 2);

	// Allocate device memory
	char *d_a, *d_b;
	// Allocate device-global score 'linear matrix'
	int *d_score;
	cudaMalloc(&d_a, sizeof(char) * a_len);
	cudaMalloc(&d_b, sizeof(char) * b_len);
	int h_score_height = min(D, h_tileHeight) + 2;
	cudaMalloc(&d_score, sizeof(int) * (h_score_height * (MAX_SEQ_SIZE+1)));
	
	// Copy hostToDevice
	cudaMemcpy(d_a, a, sizeof(char) * a_len, cudaMemcpyHostToDevice);
	cudaMemcpy(d_b, b, sizeof(char) * b_len, cudaMemcpyHostToDevice);
	// Initialise score 'linear matrix' to zeroes
	cudaMemset(d_score, 0, sizeof(int) * (h_score_height * (MAX_SEQ_SIZE+1)));
	
	// Invoke kernel
	for (int d = 0; d < D; d++) {
		int w = (d <= a_len) ? (d + 1) : (W - d);
		//printf("(d, w) = (%d, %d)\n", d, w);
		kernel<<<blocks_per_grid, threads_per_block>>>(d_a, d_b, a_len, b_len, d_score, d, w, D);
		if ((d+1) % h_tileHeight == 0) {
			nextTile<<<blocks_per_grid, threads_per_block>>>(d_score, h_score_height, b_len);
		}
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
