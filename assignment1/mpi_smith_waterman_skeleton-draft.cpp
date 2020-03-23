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

struct Point {
    int x;
    int y;
};

std::vector<Point> locate(int antidiagIdx, int a_len, int b_len);
int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i);
// /**
//  * Helper function to compute the end address of an ordinary array
//  */
// template <typename T, size_t N> const T* endOf(const T (&a)[N]) { return a+N; }

#ifdef DEBUG
void debug_print(string name, int idx, std::vector<int> vec) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (std::vector<int>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        std::cout << *i << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int idx, std::vector<Point> vec) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (std::vector<Point>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        std::cout << i->x << ',' << i->y << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int idx, int* arr, int n) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (int i = 0; i < n; i++)
        std::cout << arr[i] << ' ';
    std::cout << std::endl;
}

void test_for(int my_rank, MPI_Comm comm) {
    std::cout << "my_rank=" << my_rank << " ";
    for (int i = 0; i < 10; i++) {
        MPI_Barrier(comm);
        std::cout << i << std::endl;
        MPI_Barrier(comm);
    }
    std::cout << std::endl;
}
#endif

int smith_waterman(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len) {
    /*
     *  Please fill in your codes here.
     */
    
    // Runs in time O(a_len + b_len)
    const int num_antidiags = a_len + b_len - 1;
    // Stores max scores of all processes in each antidiag

    // Stores max score of all processes in all antidiags by far
    int max_score = 0;

    // Input data for each antidiag
    std::vector<int> lastlastAd;
    std::vector<int> lastAd;
    // std::vector<int> lastlastAd(1, 0);
    // std::vector<int> lastAd(2, 0);
    if (my_rank == 0) {
        lastlastAd.resize(1);
        std::fill(lastlastAd.begin(), lastlastAd.end(), 0);
        lastAd.resize(2);
        std::fill(lastAd.begin(), lastAd.end(), 0);
    }

    MPI_Bcast(lastlastAd.data(), lastlastAd.size(), MPI_INT, 0, comm);
    MPI_Bcast(lastAd.data(), lastAd.size(), MPI_INT, 0, comm);


#ifdef DEBUG
    int test[10] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    int localTest[5];
    if (my_rank == 0) {
        debug_print("test from my_rank=", my_rank, test, 10);
    }
    MPI_Scatter(test, 10/p, MPI_INT, localTest, 10/p, MPI_INT, 0, comm);
    debug_print("localTest from my_rank=", my_rank, localTest, 10/p);
    test_for(my_rank, comm);
#endif

    if (my_rank == 0) {
        debug_print("lastlastAd from my_rank=", my_rank, lastlastAd);
        debug_print("lastAd from my_rank=", my_rank, lastAd);
    } else {
        debug_print("lastlastAd from my_rank=", my_rank, lastlastAd);
        debug_print("lastAd from my_rank=", my_rank, lastAd);
    }

//     for (int i = 0; i < num_antidiags; i++) {

// #ifdef DEBUG
//         std::cout << "loop begins from process " << my_rank << " at step " << i << std::endl;
// #endif


//         // Blocking wait until all processes are in sync
//         // so that work partitioning can start
//         // MPI_Barrier(comm);

//         if (my_rank != 0) {
// #ifdef DEBUG
//             std::cout << "Bcasted" << std::endl;
//             debug_print("bcasted lastAd", i, lastAd);
// #endif
//         }

//         std::vector<Point> coordinates = locate(i, a_len, b_len);

//         ////////////////////
        

// #ifdef DEBUG
//         debug_print("lastlastAd", i, lastlastAd);
//         debug_print("lastAd", i, lastAd);
//         debug_print("coordinates", i, coordinates);
// #endif
//         std::vector<int> results;

//         // compute H entries along this antidiag
//         for (int k = 0; k < coordinates.size(); k++) {
//             int h = computeH_ij(
//                     lastlastAd[k], lastAd[k], lastAd[k+1],
//                     a[coordinates[k].x], b[coordinates[k].y]
//                 );
//             max_score = std::max(max_score, h);
//             results.push_back(h);
//         }

// #ifdef DEBUG
//         debug_print("results", i, results);
// #endif

//         if (my_rank == 0) {
//             // Pop boundary zero
//             lastlastAd.assign(
//                 i < a_len - 1 ? lastAd.begin() : lastAd.begin() + 1,
//                 i < b_len - 1 ? lastAd.end()   : lastAd.end() - 1);

//             // Pad boundary zero
//             lastAd.clear();
//             if (i < a_len - 1) {
//                 lastAd.push_back(0);
//             }
//             for (auto result : results) {
//                 lastAd.push_back(result);
//             }
//             if (i < b_len - 1) {
//                 lastAd.push_back(0);
//             }
//         }
//         ////////////////////

//         // All processes send their local_max_score to root process for
//         // reducing against the maximum operator, max is stored in max_scores[i]

//         //MPI_Reduce(&local_max_score, &max_scores[i], 1, MPI_INT, MPI_MAX, 0, comm);
// //        if (my_rank == 0) {
// //            // Keep current maximum out of all antidiags
// //            // max_score = max(max_score, max_scores[i]);
// //        }

//         // MPI_Barrier(comm);
//     }
#ifdef DEBUG
    MPI_Barrier(comm);
    std::cout << "Finale from process " << my_rank << std::endl;
#endif

    // The root process returns the well received max_score
    return max_score;
}   /* smith_waterman */

// Custom datatype?
// Step T: Send (i,j) to both (i+1,j) and (i,j+1)
// Step T+1: Then (i,j) receives from (i-1,j), (i,j-1), (i-1,j-1)
//            Repeat this for as many as (row_size / p) times

// For each step in antidiagonal from 1 to A + B + 1
// Divide row_size = step+1 - max(0, step - A) - max(0, step - B) into 'p' portions
// TODO helper function
// TODO Scatter 
// For each element (i, j) in H
// input: 0 or H[i-1,j], 0 or H[i,j-1], 0 or H[i-1,j-1]
// grid: H[i, j]
// output: H[i, j] to (i-1, j), (i, j-1), (i-1, j-1)


// int sizeOfAntidiag(int antidiagIdx, int a_len, int b_len) {
//     return antidiagIdx + 1 + max(0, antidiagIdx - a_len) + max(0, antidiagIdx - b_len);
// }

/**
 * Returns all the coordinates (x, y) that needs to be calculated: H[x, y] in this antidiag
 */
std::vector<Point> locate(int antidiagIdx, int a_len, int b_len) {
    std::vector<Point> coordinates;
    for (int x = 0; x <= antidiagIdx; x++) {
        int y = antidiagIdx - x;
        if (x < a_len && y < b_len) {
            Point p { x, y };
            coordinates.push_back(p);
        }
    }
    return coordinates;
}


/**
 * lastlastAd of size R,
 * lastAd of size R+1
 */
int computeH_ij(int lastlast, int last_L, int last_R, char a_i, char b_i) {
    int candidates[] = { 0, lastlast + sub_mat(a_i, b_i), last_L - GAP, last_R - GAP };
    return *std::max_element(candidates, candidates + 4);
}

// TODO scatterv allgather; locals 'task size'
/*
int* indices = NULL;
        if (my_rank == 0) {
            indices = new int[antidiagSize];
            // TODO: fill in indices to be split
            for (int j = 0; j < antidiagSize; j++) {
                indices[j] = j;
            }
            MPI_Scatter(&indices, taskSize, MPI_INT, 
                    localAntidiag, taskSize, MPI_INT, 
                    0, comm);
            delete[] indices;
        } else {
            MPI_Scatter(&indices, taskSize, MPI_INT, 
                    localAntidiag, taskSize, MPI_INT, 
                    0, comm);
        }
*/
// TODO handle p not divisible case
// TODO consider std::transform