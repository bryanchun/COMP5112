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
        a = new char[a_len];
        b = new char[b_len];
    }
    MPI_Bcast(a, a_len+1, MPI_INT, 0, comm);
    MPI_Bcast(b, b_len+1, MPI_INT, 0, comm);
    
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

// void computeHx_grid(int width, int height, int aIdx, int bIdx, int grid_max_scores[]) {

// }

// Returns local_max_score
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int local_n, int bIdx) {
    /* local-to-process variables */
#ifdef DEBUG
    std::cout << "enter computeHx" << std::endl;
#endif
    int local_max_scores[local_n];
    for (int j = 0; j < local_n; j++) {
        local_max_scores[j] = 0;
    }
#ifdef DEBUG
    std::cout << "assigned local_max_scores" << std::endl;
#endif

    // int grid_width = 1000;
    // int grid_height = 1000;
    // int num_batch_w = local_n / batch_width + (local_n % batch_width == 0 ? 0 : 1);

    int local_hx[a_len][local_n];        /* out */
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
    // for each grid h
    //      for each grid w
    //          computeHx_grid()
    //          return intercol[height], interrow[width]

    for (int i = 0; i < a_len; i++) {
        for (int j = 0; j < local_n; j++) {
            if (my_rank > 0 && j == 0) {                                                        // if not leftmost process and is leftmost local process, blocking receive
                MPI_Recv(&prev_hx[i], 1, MPI_INT, my_rank - 1,
                    my_rank - 1, comm, MPI_STATUS_IGNORE);
#ifdef DEBUG
                debug_print("prev_hx", my_rank, i, prev_hx, i+1);
#endif
            }
            local_hx[i][j] = computeH_ij(
                /* lastlast */
                i == 0 || (my_rank == 0 && j == 0)
                    ? 0 
                    : (my_rank > 0 && j == 0 
                        ? prev_hx[i-1] 
                        : local_hx[i-1][j-1]),
                    // if is top row or is first col in root process -> use 0
                    // else if first col in a non-root process -> use received col
                    // otherwise -> use last local col
                /* last_L */
                (my_rank == 0 && j == 0)
                    ? 0
                    : (my_rank > 0 && j == 0 
                        ? prev_hx[i]
                        : local_hx[i][j-1]),
                    // if is first col in root process -> use 0
                    // else if first col in a non-root process -> use received col
                    // otherwise -> use last local col
                /* last_R */
                i == 0
                    ? 0
                    : local_hx[i-1][j],
                    // if top row -> use 0
                    // else -> use last local col
                a[i], b[bIdx + j]);
#ifdef DEBUG
            debug_print2d("local_hx", my_rank, i, j, local_hx[i], j+1);
#endif
            if (my_rank < p - 1 && j == local_n - 1) {                             // if not rightmost process, send
                MPI_Send(&local_hx[i][j], 1, MPI_INT, my_rank + 1,
                    my_rank, comm);
            }

            local_max_scores[j] = max(local_max_scores[j], local_hx[i][j]);
        }
    }
    return *std::max_element(local_max_scores, local_max_scores + local_n);
}