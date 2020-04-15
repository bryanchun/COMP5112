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
void debug_print(string name, int my_rank, int h, int w, int arr[], int n) {
    std::cout << my_rank << " " << name << " in (w, h): " << w << ", " << h << " is: ";
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
// Inputs
int threads_count;
char *A, *B;
int A_len, B_len;
// Intermediates
int n_diagonal;
int** scores;
// sem_t* local_sems;
pthread_barrier_t barrier;
// sem_t count_sem;
// pthread_mutex_t counter_lock;
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
    n_diagonal = a_len + num_threads - 1;
    // local_sems = new sem_t[n_diagonal];
    // sem_init(&local_sems[0], 0, a_len);
    // for (int thread = 0; thread < n_diagonal; thread++) {
    //     sem_init(&local_sems[thread], 0, 0);
    // }
    pthread_barrier_init(&barrier, NULL, threads_count);
    // sem_init(&count_sem, 0, 1);
    // pthread_mutex_init(&counter_lock, NULL);

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
    cout << "scores:" << endl;
    for (int i = 0; i < a_len + 1; i++) {
        for (int j = 0; j < b_len + 1; j++) {
            cout << scores[i][j] << " ";
        }
        cout << endl;
    }
#endif

    // Teardown
    for (int i = 0; i <= a_len; i++) {
        delete[] scores[i];
    }
    delete[] scores;

    // for (int thread = 0; thread < n_diagonal; thread++) {
    //     sem_destroy(&local_sems[thread]);
    // }
    // delete[] local_sems;
    // sem_destroy(&count_sem);
    // pthread_mutex_destroy(&counter_lock);
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
    int base_width = B_len / threads_count;
    int outstanding = B_len % threads_count;
    int width = base_width
        + (my_rank < outstanding ? 1 : 0);
    int bIdx = my_rank < outstanding
        ? (my_rank * (base_width + 1))                                              // if my_rank is less than outstanding, multiply base_width+1 by my_rank
        : (outstanding * (base_width + 1) + (my_rank - outstanding) * base_width);  // otherwise, shift all outstanding base_width+1 then times as many more than outstanding with base_width

    // Compute my columns of the score matrix
    for (int d = 0; d < n_diagonal; d++) {
        int i = d - my_rank + 1;

        // sem_wait(&count_sem);    

        // Which threads need to work? Those with i > 0 && i <= A_len. Just normally proceed as no data dependency exists in this barrier
        // If i is in [1, A_len], this thread's dependencies are ready
        if (i > 0 && i <= A_len) {
#ifdef DEBUG
            printf("thread %ld is working at d = %d, i = %d\n", my_rank, d, i);
#endif            
            for (int w = 1; w <= width; w++) {
                int j = bIdx + w;
                scores[i][j] = computeH_ij(scores[i-1][j-1], scores[i][j-1], scores[i-1][j], A[i-1], B[j-1]);   // A and B uses 0-index
                local_maxes[my_rank] = max(local_maxes[my_rank], scores[i][j]);
            }
        }

        pthread_barrier_wait(&barrier);

//         pthread_mutex_lock(&mutex);
//         counter++;
//         if (counter == threads_count - 1) {
// #ifdef DEBUG
//             printf("thread %ld is posting to all at d = %d, i = %d\n", my_rank, d, i);
// #endif
//             // Resey counter for next barrier
//             counter = 0;
//             pthread_cond_broadcast(&cond);
//         } else {
// #ifdef DEBUG
//             printf("thread %ld is waiting at d = %d, i = %d, counter = %d\n", my_rank, d, i, counter);
// #endif
//             while (pthread_cond_wait(&cond, &mutex));
//         }
//         pthread_mutex_unlock(&mutex);


//         if (counter == threads_count - 1) {
//             // Until whichever last working thread finishes, it needs to post to all other threads on d's semaphore for the next barrier
// #ifdef DEBUG
//             printf("thread %ld is posting to all at d = %d, i = %d\n", my_rank, d, i);
// #endif
//             // pthread_mutex_lock(&counter_lock);
//             counter = 0;
//             // pthread_mutex_unlock(&counter_lock);
//             // sem_post(&count_sem);
//             for (int r = 0; r < threads_count - 1; r++) {
//                 sem_post(&local_sems[d]);
//             }
//         } else {
//             // All but one threads should wait on this barrier d, regardless of the thread having worked or not
//             // pthread_mutex_lock(&counter_lock);
//             counter++;
//             // pthread_mutex_unlock(&counter_lock);
//             // sem_post(&count_sem);
// #ifdef DEBUG
//             printf("thread %ld is waiting at d = %d, i = %d, counter = %d\n", my_rank, d, i, counter);
// #endif
//             sem_wait(&local_sems[d]);
//         }
    }

#ifdef DEBUG
    printf("Thread_max from rank %ld is %d\n", my_rank, local_maxes[my_rank]);
#endif
}

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}