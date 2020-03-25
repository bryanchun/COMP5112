/**
 * Name: CHUN Hiu Sang
 * Student id: 20421860
 * ITSC email: hschun@connect.ust.hk
*/

#include <algorithm>
#include <utility>
#ifdef DEBUG
#include <iostream>
#include <string>
#endif

#include "mpi_smith_waterman.h"

/*
 *  You can add helper functions and variables as you wish.
 */

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i);
// void computeHx_grid();
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int local_n, int bIdx);

#ifdef DEBUG
void debug_print(string name, int my_rank, int idx, int arr[], int n) {
    std::cout << name << "[" << idx << "]" << " from " << my_rank << " is ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int my_rank, int idx, vector<int> arr, int n) {
    std::cout << name << "[" << idx << "]" << " from " << my_rank << " is ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
void debug_print2d(string name, int my_rank, int primary, int idx, int arr[], int n) {
    std::cout << name << "[" << primary << "]" << "[" << idx << "]" << " from " << my_rank << " is ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
#endif

/**
 * Currently: p < b_len
 * 
 * TODO: Partition local H to moving frame (calculate H entries and toss away)
 *  - dynamic array?
 * TODO: send/recv once per process in a batch
 */
int smith_waterman(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {
    /*
     *  Please fill in your codes here.
     */

    /* Get inputs */
    MPI_Bcast(&a_len, 1, MPI_INT, 0, comm);
    MPI_Bcast(&b_len, 1, MPI_INT, 0, comm);
    if (my_rank > 0) {
        a = new char[a_len+1];
        b = new char[b_len+1];
    }
    MPI_Bcast(a, a_len+1, MPI_CHAR, 0, comm);
    MPI_Bcast(b, b_len+1, MPI_CHAR, 0, comm);
    
#ifdef DEBUG
    std::cout << "Get inputs done" << std::endl;
#endif

    /* Distribute work across processors */
    int base_n = b_len / p;
    int outstanding = b_len % p;
    int local_n = base_n
        + (my_rank < outstanding ? 1 : 0);                                  // processor does one more hx if rank is less than remainder
    int bIdx = my_rank < outstanding
        ? (my_rank * (base_n + 1))                                          // if my_rank is less than outstanding, multiply base_n+1 by my_rank
        : (outstanding * (base_n + 1) + (my_rank - outstanding) * base_n);  // otherwise, shift all outstanding base_n+1 then times as many more than outstanding with base_n
#ifdef DEBUG
    std::cout << my_rank << " has local_n " << local_n << std::endl;
    std::cout << my_rank << " has bIdx " << bIdx << std::endl;
#endif

    int local_max_score = computeHx(my_rank, p, comm, a, b, a_len, local_n, bIdx);

    /* Answer */
    int max_score;
    
    MPI_Reduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, 0, comm);
    return max_score;
}   /* smith_waterman */

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}

// void computeHx_frame(int width, int height, int aIdx, int bIdx, int frame_max_scores[]) {

// }

// Returns local_max_score
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int local_n, int bIdx) {
    /* local-to-process variables */
    int local_max_scores[local_n];
    for (int j = 0; j < local_n; j++) {
        local_max_scores[j] = 0;
    }

    // int frame_width = 1000;
    // int frame_height = 1000;
    // int num_frame_w = local_n / frame_width + (local_n % frame_width == 0 ? 0 : 1);

    int local_hx[local_n][a_len];        /* out */
    int prev_hx[a_len];                  /* in */

#ifdef DEBUG
    // int local_hx_test1[100][20000];
    // int local_hx_test1[1000][1000];

    // int local_hx_test1[1000][2000];
    std::cout << "declared local_hx_test1" << std::endl;
    // int local_hx_test2[20000][20000];
    // std::cout << "declared local_hx_test2" << std::endl;
    std::cout << "declared local_hx & prev_hx" << std::endl;
#endif

    // intercol[], interrow[][width]
    // for each frame h
    //      for each frame w
    //          computeHx_frame()
    //          return intercol[height], interrow[width]

    if (my_rank > 0) {                                                        // if not leftmost process and is leftmost local process, blocking receive
        MPI_Recv(&prev_hx, a_len, MPI_INT, my_rank - 1,
            my_rank - 1, comm, MPI_STATUS_IGNORE);
#ifdef DEBUG
        // debug_print("prev_hx", my_rank, 0, prev_hx, a_len);
        std::cout << "well received at " << my_rank << std::endl;
#endif
    }

    for (int i = 0; i < a_len; i++) {
        for (int j = 0; j < local_n; j++) {             
            /* lastlast */
            int lastlast = (i == 0 || (my_rank == 0 && j == 0))
                    ? 0 
                    : (my_rank > 0 && j == 0 
                        ? prev_hx[i-1] 
                        : local_hx[j-1][i-1]);
                    // if is top row or is first col in root process -> use 0
                    // else if first col in a non-root process -> use received col
                    // otherwise -> use last local col

            /* last_L */
            int last_L = (my_rank == 0 && j == 0)
                    ? 0
                    : (my_rank > 0 && j == 0 
                        ? prev_hx[i]
                        : local_hx[j-1][i]);
                    // if is first col in root process -> use 0
                    // else if first col in a non-root process -> use received col
                    // otherwise -> use last local col

            /* last_R */
            int last_R = (i == 0)
                    ? 0
                    : local_hx[j][i-1];
                    // if top row -> use 0
                    // else -> use last local col
            
            local_hx[j][i] = computeH_ij(lastlast, last_L, last_R, a[i], b[bIdx + j]);

            local_max_scores[j] = max(local_max_scores[j], local_hx[j][i]);
#ifdef DEBUG
            debug_print("local_hx", my_rank, 0, local_max_scores, local_n);
            std::cout << "local_hx at " << my_rank << "[" << j << "]" << "[" << i << "] " << local_hx[j][i] << std::endl;
#endif            
        }
    }
    if (my_rank < p - 1) {                                                       // if not rightmost process, send
        MPI_Send(&(local_hx[local_n - 1]), a_len, MPI_INT, my_rank + 1,
            my_rank, comm);
    }
#ifdef DEBUG
    debug_print("local_max_scores", my_rank, 0, local_max_scores, local_n);
    std::cout << "well sent at " << my_rank << std::endl;
#endif
    
    return *std::max_element(local_max_scores, local_max_scores + local_n);;
}