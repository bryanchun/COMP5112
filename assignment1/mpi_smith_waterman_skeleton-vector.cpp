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

std::vector<std::pair<int, int>> locate(int antidiagIdx, int a_len, int b_len);
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
void debug_print(string name, int idx, std::vector<std::pair<int, int>> vec) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (std::vector<std::pair<int, int>>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        std::cout << i->first << ',' << i->second << ' ';
    std::cout << std::endl;
}
void debug_print(string name, int idx, int arr[]) {
    std::cout << name << "[" << idx << "]" << ' ';
    for (int i = 0; i < sizeof(arr)/sizeof(int); i++)
        std::cout << arr[i] << ' ';
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
    std::vector<int> max_scores(num_antidiags,  0); 
    // Stores max score of all processes in all antidiags by far
    int total_max_score = 0;

    // Input data for each antidiag
    std::vector<int> lastlastAd(1, 0);
    std::vector<int> lastAd(2, 0);
    // std::vector<int> lastlastAd(1, 0);
    // std::vector<int> lastAd(2, 0);

    for (int i = 0; i < num_antidiags; i++) {
        // Blocking wait until all processes are in sync
        // so that work partitioning can start
        MPI_Barrier(comm);

        ////////////////////
        std::vector<std::pair<int, int>> coordinates = locate(i, a_len, b_len);

#ifdef DEBUG
        debug_print("lastlastAd", i, lastlastAd);
        debug_print("lastAd", i, lastAd);
        debug_print("coordinates", i, coordinates);
#endif
        std::vector<int> results;

        // compute H entries along this antidiag
        for (int k = 0; k < coordinates.size(); k++) {
            int h = computeH_ij(
                    lastlastAd[k], lastAd[k], lastAd[k+1],
                    a[coordinates[k].first], b[coordinates[k].second]
                );
            max_scores[i] = std::max(max_scores[i], h);
            results.push_back(h);
        }

#ifdef DEBUG
        debug_print("max_scores", i, max_scores);
        debug_print("results", i, results);
#endif

        // Pop boundary zero
        lastlastAd.assign(
            i < a_len - 1 ? lastAd.begin() : lastAd.begin() + 1,
            i < b_len - 1 ? lastAd.end()   : lastAd.end() - 1);

        // Pad boundary zero
        lastAd.clear();
        if (i < a_len - 1) {
            lastAd.push_back(0);
        }
        for (auto result : results) {
            lastAd.push_back(result);
        }
        if (i < b_len - 1) {
            lastAd.push_back(0);
        }

        // lastAd.assign(
        //     results.size()
        //     + i < a_len - 1 ? 1 : 0
        //     + i < b_len - 1 ? 1 : 0, 0);
        // std::copy(results.begin(), results.end(), 
        //     i < a_len - 1 ? lastAd.begin() + 1 : lastAd.begin());
        ////////////////////

        // All processes send their local_max_score to root process for
        // reducing against the maximum operator, max is stored in max_scores[i]

        //MPI_Reduce(&local_max_score, &max_scores[i], 1, MPI_INT, MPI_MAX, 0, comm);
        if (my_rank == 0) {
            // Keep current maximum out of all antidiags
            total_max_score = max(total_max_score, max_scores[i]);
        }
    }

    // The root process returns the well received total_max_score
    return total_max_score;
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
std::vector<std::pair<int, int>> locate(int antidiagIdx, int a_len, int b_len) {
    std::vector<std::pair<int, int>> coordinates;
    for (int x = 0; x <= antidiagIdx; x++) {
        int y = antidiagIdx - x;
        if (x < a_len && y < b_len) {
            coordinates.push_back(std::make_pair(x, y));
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