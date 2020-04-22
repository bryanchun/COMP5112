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
#include <sstream>
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
const int FRAME_SIZE = 10;
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
pthread_mutex_t mutex;
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

    frame_has_trailing = a_len % FRAME_SIZE > 0;
    frames_count = (a_len / FRAME_SIZE) + (frame_has_trailing ? 1 : 0);

    n_diagonal = frames_count + num_threads - 1;
    pthread_barrier_init(&barrier, NULL, threads_count);

    pthread_mutex_init(&mutex, NULL);

    // Threading
    pthread_t* thread_handles = new pthread_t[num_threads];

    for (long thread = 0; thread < num_threads; thread++) {
      pthread_create(&thread_handles[thread], NULL,
          Thread_max, (void*) thread);
    }

    for (long thread = 0; thread < num_threads; thread++) {
        pthread_join(thread_handles[thread], NULL);
    }

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

    // Get my partition
    int my_width = widths[my_rank];
    int bIdx = my_rank < outstanding
        ? (my_rank * (base_width + 1))                                              // if my_rank is less than outstanding, multiply base_width+1 by my_rank
        : (outstanding * (base_width + 1) + (my_rank - outstanding) * base_width);  // otherwise, shift all outstanding base_width+1 then times as many more than outstanding with base_width

    // Answer
    int local_max = 0;

    // Local scores is (a_len + 1) x (my_width)
    int** scores = new int*[A_len + 1];
    for (int i = 0; i <= A_len; i++) {
        scores[i] = new int[my_width]();
    }    
    
    // Compute my columns of the score matrix
    for (int d = 0; d < n_diagonal; d++) {
        int f = d - my_rank + 1;
        int my_height, aIdx;

        pthread_barrier_wait(&barrier);

        // Which threads need to work? Those with i > 0 && i <= A_len. Just normally proceed as no data dependency exists in this barrier
        // If i is in [1, A_len], this thread's dependencies are ready
        if (f > 0 && f <= frames_count) { 
            
            my_height = (f == frames_count && frame_has_trailing) ? (A_len % FRAME_SIZE) : FRAME_SIZE;
            aIdx = (f-1) * FRAME_SIZE;
            for (int h = 1; h <= my_height; h++) {
                int i = aIdx + h;               // i is 1-indexed (row 0 is padded with 0)
                for (int w = 0; w < my_width; w++) {
                    int j = bIdx + w;           // j is 0-indexed (col -1 is 0 for my_rank == 0 or buffer otherwise) - buffer is 0-indexed too
                    int lastlast = (w == 0) ? 
                                    ((my_rank == 0) ? 0 : (buffer[my_rank-1][h-1])) :      // buffer use 0-indexed, i is 1-indexed
                                    (scores[i-1][w-1]);
                    int lastL = (w == 0) ? 
                                    ((my_rank == 0) ? 0 : (buffer[my_rank-1][h])) :
                                    (scores[i][w-1]);
                    int lastR = scores[i-1][w];
                    scores[i][w] = computeH_ij(lastlast, lastL, lastR, A[i-1], B[j]);   // A and B use 0-indexed, i is 1-indexed, j is 1-indexed
#ifdef DEBUG
                    printf("thread %ld is working at d = %d, f = %d, i = %d, j = %d: lastlast = %d, lastL = %d, lastR = %d, a = %c, b = %c => score[i][w] = %d\n",
                        my_rank, d, f, i, j, lastlast, lastL, lastR, A[i-1], B[j], scores[i][w]);
#endif
                    local_max = max(local_max, scores[i][w]);
                }
            }

            if (f == frames_count) {
                local_maxes[my_rank] = local_max;
            }
        }

        pthread_barrier_wait(&barrier);

        if (f > 0 && f <= frames_count && threads_count > 1 && my_rank < threads_count - 1) {
            pthread_mutex_lock(&mutex);
            // Carry diagonal on
            buffer[my_rank][0] = scores[aIdx][my_width - 1];
            // Pass to buffer
            for (int h = 1; h <= my_height; h++) {
                int i = aIdx + h;
                buffer[my_rank][h] = scores[i][my_width - 1];
            }
            pthread_mutex_unlock(&mutex);
        }
    }

#ifdef DEBUG
    stringstream msg;
    msg << "Thread " << my_rank << " has scores:\n";
    for (int i = 0; i <= A_len; i++) {
        for (int j = 0; j < my_width; j++) {
            msg << scores[i][j] << "\t\t";
        }
        msg << endl;
    }
    cout << msg.str();
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