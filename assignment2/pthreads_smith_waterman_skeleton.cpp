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
#include <iostream>
void debug_print(string name, int my_rank, int arr[], int n) {
    std::cout << my_rank << " " << name << " is: ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int my_rank, int** arr, int n, int m) {
    std::cout << my_rank << " " << name << " is: " << endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            std::cout << arr[i][j] << ' ';
        }
        std::cout << std::endl;
    }
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
const int FRAME_SIZE = 100;
// Inputs
int threads_count;
char *A, *B;
int A_len, B_len;
int frames_count;
bool frame_has_trailing;
// Intermediates
int** scores;
sem_t* local_sems;
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

    scores = new int*[a_len+1];
    for (int i = 0; i <= a_len; i++) {
        scores[i] = new int[b_len+1]();
    }

    frame_has_trailing = a_len % FRAME_SIZE > 0;
    frames_count = (a_len / FRAME_SIZE) + (frame_has_trailing ? 1 : 0);
    local_sems = new sem_t[num_threads];
    sem_init(&local_sems[0], 0, frames_count);
    for (int thread = 1; thread < num_threads; thread++) {
        sem_init(&local_sems[thread], 0, 0);
    }

    // Threading
    pthread_t* thread_handles = new pthread_t[num_threads];

    for (long thread = 0; thread < num_threads; thread++) {
      pthread_create(&thread_handles[thread], (pthread_attr_t*) NULL,
          Thread_max, (void*) thread);
    }

    for (long thread = 0; thread < num_threads; thread++) {
        pthread_join(thread_handles[thread], NULL);
    }

#ifdef DEBUG
    // debug_print("scores", -1, scores, a_len+1, b_len+1);
#endif

    // Teardown
    for (int i = 0; i <= a_len; i++) {
        delete[] scores[i];
    }
    delete[] scores;

    for (int thread = 0; thread < num_threads; thread++) {
        sem_destroy(&local_sems[thread]);
    }
    delete[] local_sems;

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
    int base_width = B_len / threads_count;
    int outstanding = B_len % threads_count;
    int width = base_width
        + (my_rank < outstanding ? 1 : 0);
    int bIdx = my_rank < outstanding
        ? (my_rank * (base_width + 1))                                              // if my_rank is less than outstanding, multiply base_width+1 by my_rank
        : (outstanding * (base_width + 1) + (my_rank - outstanding) * base_width);  // otherwise, shift all outstanding base_width+1 then times as many more than outstanding with base_width

    // Compute my columns of the score matrix
    // by dividing the parition into frames
    for (int f = 0; f < frames_count; f++) {
        // Wait for continuation signal for my left boundary
        // as many as there are
        sem_wait(&local_sems[my_rank]);

        int frame_height = (f == frames_count-1 && frame_has_trailing) ? (A_len % FRAME_SIZE) : FRAME_SIZE;
        for (int y = f*FRAME_SIZE, i = 1+y; y < f*FRAME_SIZE + frame_height; y++, i++) {
            for (int x = 0, j = 1+bIdx; x < width; x++, j++) {
                scores[i][j] = computeH_ij(scores[i-1][j-1], scores[i][j-1], scores[i-1][j], A[y], B[bIdx + x]);
                local_maxes[my_rank] = max(local_maxes[my_rank], scores[i][j]);
            }
        }

        if (my_rank < threads_count-1) {
            // Give continuation signal for right-neighbouring thread
            // as many as there are
            sem_post(&local_sems[my_rank+1]);
        }
    }

#ifdef DEBUG
    // cout << my_rank << " has called sem_wait/post _ times: " << num_calls << endl;
    cout << my_rank << " has local_max is " << local_maxes[my_rank] << endl;
#endif
}

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}