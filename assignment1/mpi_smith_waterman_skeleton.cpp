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
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len);

#ifdef DEBUG
void debug_print(string name, int idx, std::vector<int> vec) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (std::vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        std::cout << *i << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int my_rank, int idx, int* arr, int n) {
    std::cout << name << "[" << idx << "]" << " from " << my_rank << " is ";
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}
#endif

/**
 * Currently: p == b_len
 * 
 * TODO: GetInput - Bcast inputs to every process first
 * TODO: p < b_len, Scatter string b in characters
 * TODO: max(a_len, b_len)
 */
int smith_waterman(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {
    /*
     *  Please fill in your codes here.
     */

    MPI_Bcast(&a_len, 1, MPI_INT, 0, comm);
    MPI_Bcast(&b_len, 1, MPI_INT, 0, comm);
    if (my_rank > 0) {
        a = new char[a_len];
        b = new char[b_len];
    }
    for (int i = 0; i < a_len; i++) {
        MPI_Bcast(&a[i], 1, MPI_INT, 0, comm);
    }
    for (int i = 0; i < b_len; i++) {
        MPI_Bcast(&b[i], 1, MPI_INT, 0, comm);
    }
    
    /* answer */
    int max_score;
    int local_max_score = computeHx(my_rank, p, comm, a, b, a_len, b_len);

    MPI_Reduce(&local_max_score, &max_score, 1, MPI_INT, MPI_MAX, 0, comm);
    return max_score;
}   /* smith_waterman */

int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}

// Returns local_max_score
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {
    /* local-to-process variables */
    int local_max_score = 0;    /* out */
    int local_hx[a_len];        /* out */
    int prev_hx[a_len];         /* in */

    for (int i = 0; i < a_len; i++) {
        if (my_rank > 0) {                                  // if not leftmost progress, blocking receive
            MPI_Recv(&prev_hx[i], 1, MPI_INT, my_rank - 1,
                my_rank - 1, comm, MPI_STATUS_IGNORE);
#ifdef DEBUG
            debug_print("prev_hx", my_rank, i, prev_hx, i+1);
#endif
        }
        local_hx[i] = computeH_ij(
            my_rank > 0 && i > 0    ? prev_hx[i-1]  : 0,    // if leftmost process or top row, use 0
            my_rank > 0             ? prev_hx[i]    : 0,    // if leftmost process, use 0
            i > 0                   ? local_hx[i-1] : 0,    // if top row, use 0
            a[i], b[my_rank]);
#ifdef DEBUG
        debug_print("local_hx", my_rank, i, local_hx, i+1);
#endif
        if (my_rank < p - 1) {                             // if not rightmost process, send
            MPI_Send(&local_hx[i], 1, MPI_INT, my_rank + 1,
                my_rank, comm);
        }

        local_max_score = max(local_max_score, local_hx[i]);
    }
    return local_max_score;
}