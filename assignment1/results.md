| p       | Input     | Hostfile | Ans  | Running time |
| ------- | --------- | -------- | ---- | ------------ |
| serial  | 2k.in     |          |      |              |
| **1**   | 2k.in     |          |      |              |
| **2**   | 2k.in     |          |      |              |
| serial  | sample.in |          | 13   | 2.1828e-05 s |
| pre-mpi | sample.in |          | 13   | 4.1813e-05 s |
| **4**   | sample.in | Yes      |      |              |
| serial  | 1k.in     |          |      |              |
| **4**   | 1k.in     | Yes      |      |              |
| **7**   | 2k.in     | Yes      |      |              |



4. ```
	csl2wk34:hschun:171> ./serial_smith_waterman ../datasets/sample.in
           T       G       T       T       A       C       G   G
   G       0       3       1       0       0       0       3   3
   G       0       3       1       0       0       0       3   6
   T       3       1       6       4       2       0       1   4
   T       3       1       4       9       7       5       3   2
   G       1       6       4       7       6       4       8   6
   A       0       4       3       5       10      8       6   5
   C       0       2       1       3       8       13      11  9
   T       3       1       5       4       6       11      10  8
   A       1       0       3       2       7       9       8   7
   13
   Time: 0.00112812 s
  ```
  ```
  csl2wk34:hschun:162> mpiexec -n 1 ./mpi_smith_waterman datasets/sample.in
  lastlastAd[0] 0 
  lastAd[0] 0 0 
  coordinates[0] 0,0 
  max_scores[0] 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
  results[0] 0 
  lastlastAd[1] 0 0 
  lastAd[1] 0 0 0 
  coordinates[1] 0,1 1,0 
  max_scores[1] 0 3 0 0 0 0 0 0 0 0 0 0 0 0 0 0 
  results[1] 3 0 
  lastlastAd[2] 0 0 0 
  lastAd[2] 0 3 0 0 
  coordinates[2] 0,2 1,1 2,0 
  max_scores[2] 0 3 3 0 0 0 0 0 0 0 0 0 0 0 0 0 
  results[2] 1 3 3 
  lastlastAd[3] 0 3 0 0 
  lastAd[3] 0 1 3 3 0 
  coordinates[3] 0,3 1,2 2,1 3,0 
  max_scores[3] 0 3 3 3 0 0 0 0 0 0 0 0 0 0 0 0 
  //results[3] 0 1 1 3 
  lastlastAd[4] 0 1 3 3 0 
  lastAd[4] 0 0 1 1 3 0 
  coordinates[4] 0,4 1,3 2,2 3,1 4,0 
  max_scores[4] 0 3 3 3 6 0 0 0 0 0 0 0 0 0 0 0 
  results[4] 0 0 6 1 1 
  lastlastAd[5] 0 0 1 1 3 0 
  lastAd[5] 0 0 0 6 1 1 0 
  coordinates[5] 0,5 1,4 2,3 3,2 4,1 5,0 
  max_scores[5] 0 3 3 3 6 6 0 0 0 0 0 0 0 0 0 0 
  results[5] 0 0 4 4 6 0 
  lastlastAd[6] 0 0 0 6 1 1 0 
  lastAd[6] 0 0 0 4 4 6 0 0 
  coordinates[6] 0,6 1,5 2,4 3,3 4,2 5,1 6,0 
  max_scores[6] 0 3 3 3 6 6 9 0 0 0 0 0 0 0 0 0 
  results[6] 3 0 2 9 4 4 0 
  lastlastAd[7] 0 0 0 4 4 6 0 0 
  lastAd[7] 0 3 0 2 9 4 4 0 0 
  coordinates[7] 0,7 1,6 2,5 3,4 4,3 5,2 6,1 7,0 
  max_scores[7] 0 3 3 3 6 6 9 7 0 0 0 0 0 0 0 0 
  results[7] 3 3 0 7 7 3 2 3 
  lastlastAd[8] 0 3 0 2 9 4 4 0 0 // 3 0 2 9 4 4 0 0 
  lastAd[8] 0 3 3 0 7 7 3 2 3 	// 3 3 0 7 7 3 2 3 0
  coordinates[8] 1,7 2,6 3,5 4,4 5,3 6,2 7,1 8,0 
  max_scores[8] 0 3 3 3 6 6 9 7 6 0 0 0 0 0 0 0 
  results[8] 3 1 1 5 6 5 1 1 						// 6 1 5 6 5 1 1 1
  lastlastAd[9] 0 3 3 0 7 7 3 2 3 
  lastAd[9] 3 1 1 5 6 5 1 1 
  coordinates[9] 2,7 3,6 4,5 5,4 6,3 7,2 8,1 
  max_scores[9] 0 3 3 3 6 6 9 7 6 10 0 0 0 0 0 0 
  results[9] 1 0 3 4 4 10 0 						// 4 3 10 5 3 5 0
  lastlastAd[10] 3 1 1 5 6 5 1 1 
  lastAd[10] 1 0 3 4 4 10 0 
  coordinates[10] 3,7 4,6 5,5 6,4 7,3 8,2 
  max_scores[10] 0 3 3 3 6 6 9 7 6 10 9 0 0 0 0 0 
  results[10] 0 4 2 2 9 8  						  // 2 8 8 8 4 3
  lastlastAd[11] 1 0 3 4 4 10 0 
  lastAd[11] 0 4 2 2 9 8 
  coordinates[11] 4,7 5,6 6,5 7,4 8,3 
  max_scores[11] 0 3 3 3 6 6 9 7 6 10 9 7 0 0 0 0 
  results[11] 4 2 6 7 7 								// 6 6 13 6 2
  lastlastAd[12] 0 4 2 2 9 8 
  lastAd[12] 4 2 6 7 7 
  coordinates[12] 5,7 6,6 7,5 8,4 
  max_scores[12] 0 3 3 3 6 6 9 7 6 10 9 7 5 0 0 0 
  results[12] 2 4 5 5 									// 5 11 11 7
  lastlastAd[13] 4 2 6 7 7 
  lastAd[13] 2 4 5 5 
  coordinates[13] 6,7 7,6 8,5 
  max_scores[13] 0 3 3 3 6 6 9 7 6 10 9 7 5 3 0 0 
  results[13] 2 3 3 										// 9 10 9
  lastlastAd[14] 2 4 5 5 
  lastAd[14] 2 3 3 
  coordinates[14] 7,7 8,6 
  max_scores[14] 0 3 3 3 6 6 9 7 6 10 9 7 5 3 1 0 
  results[14] 1 1 											// 8 8
  lastlastAd[15] 2 3 3 
  lastAd[15] 1 1 
  coordinates[15] 8,7 
  max_scores[15] 0 3 3 3 6 6 9 7 6 10 9 7 5 3 1 0 
  results[15] 0 												// 7
  10
  Time: 0.000691229 s
  ```
  
  