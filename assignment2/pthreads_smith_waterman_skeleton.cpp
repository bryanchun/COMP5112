/**
 * Name: CHUN Hiu Sang
 * Student id: 20421860
 * ITSC email: hschun@connect.ust.hk
*/

#include <pthread.h>
#include <semaphore.h>
#include <algorithm>
#include "pthreads_smith_waterman.h"
#ifdef DEBUG
#include <stdio.h>
#include <iostream>
void debug_print(string name, int my_rank, int arr[], int n) {
    std::cout << my_rank << " " << name << " is: ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
#endif

/*
 *  You can add helper functions and variables as you wish.
 */
void* Thread_max(void* rank);
int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i);

/* Global variables for communication */
// Constants
const int MAX_THREADS = 8;
const int FRAME_SIZE = 1;
// Inputs
int threads_count;
char *A, *B;
int A_len, B_len;
int frames_count;
bool frame_has_trailing;
int base_width;
int outstanding;
int* widths;
// Intermediates
int n_diagonal;
int** buffer;
// int*** scores;
pthread_barrier_t barrier;
// Outputs
int local_maxes[MAX_THREADS] = {};

int smith_waterman(int num_threads, char *a, char *b, int a_len, int b_len){
    /*
     *  Please fill in your codes here.
     */
    
    // Initialise global variables
    threads_count = num_threads;
    A = a; B = b;
    A_len = a_len; B_len = b_len;
    
    base_width = b_len / threads_count;
    outstanding = b_len % threads_count;

    widths = new int[threads_count];
    for (int r = 0; r < threads_count; r++) {
        widths[r] = base_width + (r < outstanding ? 1 : 0);
    }

    // Producer-Consumer buffer for more threads
    if (threads_count > 1) {
        buffer = new int*[threads_count - 1];
        for (int r = 0; r < threads_count - 1; r++) {
            buffer[r] = new int[FRAME_SIZE + 1]();
        }
    }

#ifdef DEBUG
    cout << "widths:" << endl;
    for (int r = 0; r < threads_count; r++) {
        cout << widths[r] << " ";
    }
    cout << endl;
    // cout << "scores:" << endl;
    // for (int i = 0; i <= a_len; i++) {
    //     for (int r = 0; r < threads_count; r++) {
    //         for (int j = 0; j < widths[r]; j++) {
    //             cout << scores[r][i][j] << " ";
    //         }
    //     }
    //     cout << endl;
    // }
#endif

    frame_has_trailing = a_len % FRAME_SIZE > 0;
    frames_count = (a_len / FRAME_SIZE) + (frame_has_trailing ? 1 : 0);

    n_diagonal = frames_count + num_threads - 1;
    pthread_barrier_init(&barrier, NULL, threads_count);

    // Threading
    pthread_t* thread_handles = new pthread_t[num_threads];

    for (long thread = 0; thread < num_threads; thread++) {
      pthread_create(&thread_handles[thread], NULL,
          Thread_max, (void*) thread);
    }

    for (long thread = 0; thread < num_threads; thread++) {
        pthread_join(thread_handles[thread], NULL);
    }

#ifdef DEBUG
    // cout << "scores:" << endl;
    // for (int i = 0; i <= a_len; i++) {
    //     for (int r = 0; r < threads_count; r++) {
    //         for (int j = 0; j < widths[r]; j++) {
    //             cout << scores[r][i][j] << " ";
    //         }
    //     }
    //     cout << endl;
    // }
#endif

    // Teardown

    if (threads_count > 1) {
        for (int r = 0; r < threads_count - 1; r++) {
            delete[] buffer[r];
        }
        delete[] buffer;
    }

    delete[] widths;
    
    pthread_barrier_destroy(&barrier);

    delete[] thread_handles;

    // Reduce and return answer
#ifdef DEBUG
    debug_print("local_maxes", -1, local_maxes, num_threads);
#endif
	return *std::max_element(local_maxes, local_maxes + num_threads);
}

void* Thread_max(void* rank) {
    long my_rank = (long) rank;
    int my_width = widths[my_rank];

    // Local scores is (a_len + 1) x (my_width + 1)
    int** scores = new int*[A_len + 1];
    for (int i = 0; i <= A_len; i++) {
        scores[i] = new int[my_width + 1]();
    }

    // Get my partition
    int bIdx = my_rank < outstanding
        ? (my_rank * (base_width + 1))                                              // if my_rank is less than outstanding, multiply base_width+1 by my_rank
        : (outstanding * (base_width + 1) + (my_rank - outstanding) * base_width);  // otherwise, shift all outstanding base_width+1 then times as many more than outstanding with base_width

    // Compute my columns of the score matrix
    for (int d = 0; d < n_diagonal; d++) {
        int f = d - my_rank + 1;

        // TODO: Initialize leftmost col by getting from last

        // Which threads need to work? Those with i > 0 && i <= A_len. Just normally proceed as no data dependency exists in this barrier
        // If i is in [1, A_len], this thread's dependencies are ready
        if (f > 0 && f <= frames_count) {
#ifdef DEBUG
            // printf("thread %ld is working at d = %d, f = %d\n", my_rank, d, f);
#endif            
            int frame_height = (f == frames_count && frame_has_trailing) ? (A_len % FRAME_SIZE) : FRAME_SIZE;
            int aIdx = 1 + (f-1) * FRAME_SIZE;
            for (int b = 0, i = aIdx; b < frame_height, i < aIdx + frame_height; b++, i++) {
                for (int w = 0; w < my_width; w++) {
                    int j = 1 + bIdx + w;
                    int lastlast = (j == 1) ? 
                                    ((my_rank == 0) ? 0 : buffer[my_rank-1][b-2]):      // buffer use 0-indexed, i is 1-indexed
                                    (scores[i-1][j-1]);
                    int lastL = (j == 1) ? 
                                    ((my_rank == 0) ? 0 : buffer[my_rank-1][b-1]):
                                    (scores[i][j-1]);
                    int lastR = scores[i-1][j];
#ifdef DEBUG
                    printf("thread %ld is working at d = %d, f = %d, i = %d, j = %d: lastlast = %d, lastL = %d, lastR = %d\n",
                        my_rank, d, f, i, j, lastlast, lastL, lastR);
#endif
                    scores[i][j] = computeH_ij(lastlast, lastL, lastR, A[i-1], B[j-1]);   // A and B use 0-indexed, i is 1-indexed, j is 1-indexed
                    local_maxes[my_rank] = max(local_maxes[my_rank], scores[i][j]);
                }

                // Pass to buffer
                buffer[my_rank][b] = scores[i][my_width - 1];
            }
        }

        pthread_barrier_wait(&barrier);
#ifdef DEBUG
        if (my_rank < threads_count - 1) {
            printf("Thread %ld has buffer:\n", my_rank);
            for (int i = 0; i < FRAME_SIZE + 1; i++) {
                cout << buffer[my_rank][i] << " ";
            }
            cout << endl;
        }
#endif
    }

#ifdef DEBUG
    printf("Thread %ld has scores:\n", my_rank);
    for (int i = 0; i <= A_len; i++) {
        for (int j = 0; j <= my_width; j++) {
            cout << scores[i][j] << " ";
        }
        cout << endl;
    }
#endif

    // Teardown
    for (int i = 0; i <= A_len; i++) {
        delete[] scores[i];
    }
    delete[] scores;

#ifdef DEBUG
    printf("Thread_max from rank %ld is %d\n", my_rank, local_maxes[my_rank]);
#endif
}

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}