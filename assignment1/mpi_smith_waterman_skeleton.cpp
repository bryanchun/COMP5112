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
vector<int> computeHx_frame(
    int left_right_frame_buffer[], int top_down_frame_buffer[], int* lower_right_to_upper_left_buffer,
    int width, int height, char* a, char* b, int aIdx, int bIdx);
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int local_n, int bIdx);

#ifdef DEBUG
void debug_print(string name, int my_rank, int idx, int arr[], int n) {
    std::cout << my_rank << " has " << name << "[" << idx << "]" << " is ";
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
 * TODO: truncate at non-divisible frames
 * TODO: handle small inputs still
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

// computeHx_frame(left_right_frame_buffer, top_down_frame_buffer[w-1], lower_right_to_upper_left_buffer)
// Update left_right_frame_buffer pass by ref
// Update top_down_frame_buffer pass by ref
// Update lower_right_to_upper_left_buffer pass by ref
// Return frame_max_scores[width] as a vector
vector<int> computeHx_frame(
    int left_right_frame_buffer[],           /* in/out [height] */
    int top_down_frame_buffer[],             /* in/out [width] */
    int* lower_right_to_upper_left_buffer,   /* in/out int */
    int width, int height, char* a, char* b, int aIdx, int bIdx) {


    int frame_hx[width][height];
    vector<int> frame_max_scores(width, 0);

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
#ifdef DEBUG
    cout << "computeHx_frame: setup " << i << ", " << j << endl;
#endif
        // width=1, height=2
            /* lastlast */
            int lastlast = (i == 0 && j == 0)
                ? *lower_right_to_upper_left_buffer
                : (i == 0) ? top_down_frame_buffer[j-1]
                : (j == 0) ? left_right_frame_buffer[i-1]
                : frame_hx[j-1][i-1];
                // is top-left corner -> use lower_right_to_upper_left_buffer
                // is top row -> use top_down_frame_buffer
                // is leftmost col -> use left_right_frame_buffer
                // else -> use upperleft neighbour from frame_hx
#ifdef DEBUG
    cout << "computeHx_frame: lastlast" << endl;
    cout << "left_right_frame_buffer[i] " << left_right_frame_buffer[i] << endl;
#endif
            /* last_L */
            int last_L = (j == 0)
                ? left_right_frame_buffer[i]
                : frame_hx[j-1][i];
                // is first col -> use left_right_frame_buffer
                // else -> use left neighbour from frame_hx
#ifdef DEBUG
    cout << "computeHx_frame: last_L" << endl;
#endif
            /* last_R */
            int last_R = (i == 0)
                ? top_down_frame_buffer[j]
                : frame_hx[j][i-1];
                // is first row -> use top_down_frame_buffer
                // else -> use top neighbour from frame_hx
#ifdef DEBUG
    cout << "computeHx_frame: last_R" << endl;
#endif
            frame_hx[j][i] = computeH_ij(lastlast, last_L, last_R, a[aIdx + i], b[bIdx + j]);
            frame_max_scores[j] = max(frame_max_scores[j], frame_hx[j][i]);

#ifdef DEBUG
            cout << "lastlast, last_L, last_R, a', b': " << lastlast << " " << last_L << " " << last_R << " " << a[aIdx + i] << " " << b[bIdx + j] << endl;
            cout << "frame_hx[j][i]=" << frame_hx[j][i] << endl;
            cout << "frame_max_scores[j]=" << frame_max_scores[j] << endl;
#endif 
        }
    }
#ifdef DEBUG
    cout << "frame_hx: ";
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            cout << frame_hx[j][i] << " ";
        }
        cout << endl;
    }
#endif
    for (int i = 0; i < height; i++) {
        left_right_frame_buffer[i] = frame_hx[width-1][i];
    }
    for (int j = 0; j < width; j++) {
        top_down_frame_buffer[j] = frame_hx[j][height-1];
    }
    *lower_right_to_upper_left_buffer = frame_hx[width][height];

    return frame_max_scores;
}

// Returns local_max_score
int computeHx(int my_rank, int p, MPI_Comm comm, char *a, char *b, int a_len, int b_len, int bIdx) {

    int frame_h = 1000;
    int frame_w = 1000;
    
    int num_y_full_frame = a_len / frame_h;
    int outstanding_y = a_len % frame_h;
    int has_trailing_y_frame = (outstanding_y > 0);
    int num_y_frame = outstanding_y + has_trailing_y_frame ? 1 : 0;

    int num_x_full_frame = b_len / frame_w;
    int outstanding_x = b_len % frame_w;
    int has_trailing_x_frame = (outstanding_x > 0);
    int num_x_frame = outstanding_x + has_trailing_x_frame ? 1 : 0;

    int local_max_scores[num_x_frame];
    for (int j = 0; j < num_x_frame; j++) {
        local_max_scores[j] = 0;
    }

    // Each of these store the result for this frame, and the next frame shall use it
    /* each row of frames passes 'num_x_frame' many 'frame_w' wide h entries to all frames downwards, where last frame could be less than 'frame_w' */
    int top_down_frame_buffer[num_x_frame][frame_w];
    // int** top_down_frame_buffer = new int*[num_x_frame];
    // for (int k = 0; k < num_x_frame; k++) {
    //     top_down_frame_buffer[k] = new int[frame_w];
    //     // for (int k2 = 0; k2 < frame_w; k2++) {
    //     //     top_down_frame_buffer[k][k2] = 0;
    //     // }
    // }
    /* each row of frames passes 'num_x_frame-1' many lower right corners to its right-down frame; rightmost frame of the row just send to next process later so do not bother */
    int lower_right_to_upper_left_buffer[num_x_frame];
    // int* lower_right_to_upper_left_buffer = new int[num_x_frame];
    /* each frame passes 'frame_h' high of h entries to the right frame */
    int left_right_frame_buffer[frame_h];
    // int* left_right_frame_buffer = new int[frame_h];

    for (int h = 0; h < num_y_frame; h++) {
        for (int w = 0; w < num_x_frame; w++) {

#ifdef DEBUG
    std::cout << my_rank << " is in (w, h): " << w << ", " << h << std::endl;
#endif

            int width = (w == num_x_frame - 1 && has_trailing_x_frame) ? outstanding_x : frame_w;
            int height = (h == num_y_frame - 1 && has_trailing_y_frame) ? outstanding_y : frame_h;

            /* leftmost frame of a row */
            /* Getting buffer from the left */
            if (w == 0) {
                if (my_rank > 0) {
                    // Recv left_right_frame_buffer from rank-1
                    int prevHeight;
                    MPI_Recv(&prevHeight, 1, MPI_INT, my_rank - 1,
                        my_rank - 1, comm, MPI_STATUS_IGNORE);
                    MPI_Recv(&left_right_frame_buffer, prevHeight, MPI_INT, my_rank - 1,
                        my_rank - 1, comm, MPI_STATUS_IGNORE);
                } else {
                    // Default left_right_frame_buffer to zeros
                    for (int k = 0; k < height; k++) {
                        left_right_frame_buffer[k] = 0;
                    }
                }
            }
            // else: Use left_right_frame_buffer from last frame

#ifdef DEBUG
    std::cout << my_rank << " has lrfb ready in (w, h): " << w << ", " << h << std::endl;
#endif

            /* Getting buffer from the top */
            if (h == 0) {
                // Default top_down_frame_buffer to 0
                for (int k = 0; k < width; k++) {
                    top_down_frame_buffer[w][k] = 0;
                }
            }
            // else: Use top_down_frame_buffer from last frame

#ifdef DEBUG
    std::cout << my_rank << " has tpfb ready in (w, h): " << w << ", " << h << std::endl;
#endif

            /* Getting buffer from top-left */
            if (w == 0 || h == 0) {
                // Default lower_right_to_upper_left_buffer to 0
                lower_right_to_upper_left_buffer[w] = 0;
            } 
            // else: Use lower_right_to_upper_left_buffer from last frame

#ifdef DEBUG
            std::cout << my_rank << " has lrtulb ready in (w, h): " << w << ", " << h << std::endl;
            debug_print("left_right_frame_buffer:before", h, w, left_right_frame_buffer, height);
#endif
            vector<int> frame_max_scores = computeHx_frame(
                left_right_frame_buffer,
                top_down_frame_buffer[w],
                &(lower_right_to_upper_left_buffer[w]),
                width, height, a, b,
                w * frame_w, bIdx + h * frame_h        // string indexes can only be 0, frame_w, 2*frame_w, ... etc. works for trailing frame too
            );
            local_max_scores[w] = max(local_max_scores[w], *std::max_element(frame_max_scores.begin(), frame_max_scores.end()));

#ifdef DEBUG
            debug_print("left_right_frame_buffer:after", h, w, left_right_frame_buffer, height);
            cout << "sending of height " << height << endl;
#endif

#ifdef DEBUG
            std::cout << my_rank << " has local_max_scores[w] " << local_max_scores[w] << " in (w, h): " << w << ", " << h << std::endl;
#endif

            /* rightmost frame of a row */
            /* Sending buffer to next process if not the last process */
            if (my_rank < p - 1 && w == num_x_frame - 1) {
                // Send left_right_frame_buffer to my_rank+1
                MPI_Send(&height, 1, MPI_INT, my_rank + 1,
                    my_rank, comm);
                MPI_Send(&left_right_frame_buffer, height, MPI_INT, my_rank + 1,
                    my_rank, comm);
            }
        }
        // Completed updating top_down_frame_buffer for a row
    }

    /* Tear down */
    // for (int k = 0; k < num_x_frame; k++) {
    //     delete[] top_down_frame_buffer[k];
    // }
    // delete[] top_down_frame_buffer;
    // delete[] lower_right_to_upper_left_buffer;
    // delete[] left_right_frame_buffer;

    return *std::max_element(local_max_scores, local_max_scores + num_x_frame);
}