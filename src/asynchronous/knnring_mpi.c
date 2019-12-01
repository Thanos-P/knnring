/**
 * @file   knnring_mpi.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 * This is the ASYNCRONOUS implementation of MPI
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <mpi.h>
#include "knnring.h"

// =================
// === UTILITIES ===
// =================

//! Swap the values of element a and b
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

//! Partition array elements on the respective side of the array[pivotIndex]
int partition(double *array, int *idx, int left, int right, int pivotIndex){

  double pivotValue = array[pivotIndex];
  double temp;
  SWAP(array[pivotIndex], array[right]);  // move pivot to end
  SWAP(idx[pivotIndex], idx[right]);
  int storeIndex = left;
  int i;
  for(i = left; i < right; i++){
    if(array[i] < pivotValue){
      SWAP(array[storeIndex], array[i]);
      SWAP(idx[storeIndex], idx[i]);
      storeIndex++;
    }
  }
  SWAP(array[right], array[storeIndex]);  // move pivot to its final place
  SWAP(idx[right], idx[storeIndex]);
  return storeIndex;
}

//! Put the k-th smallest element of array in index-position k
//  within left..right inclusive
void quickselect(double *array, int *idx, int left, int right, int k){

  if(left == right){      // If the list contains only one element,
    return;
  }
  // Select a random pivotIndex between left and right
  int pivotIndex = left + rand() % (right-left+1);
  pivotIndex = partition(array, idx, left, right, pivotIndex);
  // The pivot is in its final sorted position
  if(k == pivotIndex){
    return;
  } else if(k < pivotIndex){
    return quickselect(array, idx, left, pivotIndex-1, k);
  } else{
    return quickselect(array, idx, pivotIndex+1, right, k);
  }
}

// ==================
// === SEQUENTIAL ===
// ==================

//! Compute k nearest neighbors of each point in X [n-by-d]
/*!

  \param X      Corpus data points            [n-by-d]
  \param Y      Query data points             [m-by-d]
  \param n      Number of data points         [scalar]
  \param m      Number of query points        [scalar]
  \param d      Number of dimensions          [scalar]
  \param k      Number of neighbors           [scalar]

  \return The kNN result
*/
knnresult kNN(double * X, double * Y, int n, int m, int d, int k) {

  int i, j;
  // Initialize distances and indexes arrays
  double *dist = (double *)malloc(m*n*sizeof(double));
  if(dist == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }
  int *idx = (int *)malloc(m*n*sizeof(int));
  if(idx == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }

  // Find distances
  int inc = 1;
  for(i = 0; i < m; i++){
    // MATLAB: sum(Y.^2,2)
    dist[i*n] = pow(cblas_dnrm2(d, (Y+i*d), inc), 2);
    for(j = 0; j < n; j++){
      // Value to be expanded
      double row_expand;
      if(j == 0)
        row_expand = dist[i*n];
      // Expand value for all n and add MATLAB: sum(X.^2,2).'
      dist[i*n + j] = row_expand + pow(cblas_dnrm2(d, (X+j*d), inc), 2);
      idx[i*n + j] = j;
    }
  }

  enum CBLAS_ORDER order = CblasRowMajor;
  enum CBLAS_TRANSPOSE transX = CblasTrans;
  enum CBLAS_TRANSPOSE transY = CblasNoTrans;
  double alpha = -2.0, beta = 1.0;
  int ldX = d, ldY = d, lddist = n;

  // MATLAB: - 2 * Y*X.' + sum(Y.^2,2) + sum(X.^2,2).'
  cblas_dgemm(order, transY, transX, m, n, d, alpha, Y, ldY, X, ldX, beta, dist, lddist);

  // If corpus set is equal to query set, set diagonal to zero
  // (which is distance of the same point)
  // to avoid floating point errors
  if(X == Y){
    for(i = 0; i < n; i++){
      dist[i*n + i] = 0.0;
    }
  }

  // MATLAB: sqrt(ans)
  for(i = 0; i < n*m; i++){
    dist[i] = sqrt(dist[i]);
  }

  // Initialize result
  knnresult result;
  result.m = m;
  result.k = k;
  result.nidx = (int *)malloc(m*k*sizeof(int));
  if(result.nidx == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }
  result.ndist = (double *)malloc(m*k*sizeof(double));
  if(result.ndist == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }

  // Find k nearest neighbors using quickselect for each query point
  for(i = 0; i < m; i++){
    for(j = 0; j < k; j++){
      quickselect((dist+i*n), (idx+i*n), 0, n-1, j);
      // Store them in the appropriate positions of the struct
      result.nidx[i*k + j] = idx[i*n + j];
      result.ndist[i*k + j] = dist[i*n + j];
    }
  }

  free(dist);
  free(idx);

  return result;
}

// ===========
// === MPI ===
// ===========

//! Merge kNN for new corpus points with the kNN found previously
/*!

  \param selfknn        old kNN to be merged
  \param newknn         kNN found in last search

*/
void mergekNN(knnresult *selfknn, knnresult *newknn){

  knnresult result;
  int selfcount, newcount, i, j;
  // Initialize merged kNN result
  result.m = selfknn->m;
  result.k = selfknn->k;
  result.nidx = (int *)malloc(result.m*result.k*sizeof(int));
  if(result.nidx == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }
  result.ndist = (double *)malloc(result.m*result.k*sizeof(double));
  if(result.ndist == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }

  // For every query point
  for(i = 0; i < result.m; i++){
    selfcount = 0;
    newcount = 0;
    // Compare nearest neighbor of each search and add the nearest
    for(j = 0; j < result.k; j++){
      if(selfknn->ndist[i*result.k + selfcount] <= newknn->ndist[i*result.k + newcount]){
        result.ndist[i*result.k + j] = selfknn->ndist[i*result.k + selfcount];
        result.nidx[i*result.k + j] = selfknn->nidx[i*result.k + selfcount];
        selfcount++;
      }else{
        result.ndist[i*result.k + j] = newknn->ndist[i*result.k + newcount];
        result.nidx[i*result.k + j] = newknn->nidx[i*result.k + newcount];
        newcount++;
      }
    }
  }
  // Update selfknn to be the merged result
  free(selfknn->ndist);
  free(selfknn->nidx);
  selfknn->ndist = result.ndist;
  selfknn->nidx = result.nidx;
}

//! Compute distributed all-kNN of points in X
/*!

  \param X      Data points                [n-by-d]
  \param n      Number of data points      [scalar]
  \param d      Number of dimensions       [scalar]
  \param k      Number of neighbors        [scalar]

  \return The kNN result
*/
knnresult distrAllkNN(double * X, int n, int d, int k){

  int tid, numTasks;
  // Error handling variables
  int err, errlen;
  char errbuffer[MPI_MAX_ERROR_STRING];
  int i, j;
  MPI_Status *mpistat = (MPI_Status *)malloc(2 * sizeof(MPI_Status));
  MPI_Request *mpireq = (MPI_Request*)malloc(2 * sizeof(MPI_Request));

  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &tid);

  if(tid == 0){
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
    printf("Running asynchronous implementation\n");
    printf("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
  }

  // Data arrays / structs
  double *corpusbuffer = (double *)malloc(n*d*sizeof(double));
  if(corpusbuffer == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }
  double *sendbuffer = (double *)malloc(n*d*sizeof(double));
  if(sendbuffer == NULL){
    printf("Failed to allocate memory\n");
    exit(EXIT_FAILURE);
  }

  // selfknn is the result that is returned
  // newknn tracks new neighbors each iteration
  knnresult selfknn, newknn;

  // Send / receive info
  int dst = (tid+1)%numTasks;
  int rcv = (tid-1+numTasks)%numTasks;
  int tag = 1;
  int size = n*d;

  // Hide communication costs by sending/recieving while computing
  // Send query points of each process to next process
  err = MPI_Isend(X, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, mpireq);
  if(err){
    MPI_Error_string(err, errbuffer, &errlen);
    printf("Error %d [%s]\n", err, errbuffer);
  }
  // Receive new points as corpus points and keep query points for new search
  err = MPI_Irecv(corpusbuffer, size, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, mpireq+1);
  if(err){
    MPI_Error_string(err, errbuffer, &errlen);
    printf("Error %d [%s]\n", err, errbuffer);
  }

  // First iteration for the points of each process
  selfknn = kNN(X, X, n, n, d, k);
  // Offset indexes to match starting array
  for(i = 0; i < selfknn.m * selfknn.k; i++){
    selfknn.nidx[i] += (tid-1+numTasks)%numTasks * n;
  }

  // Wait for sends/receives to finish before sending/recieving again
  MPI_Waitall(2, mpireq, mpistat);

  for(i = 0; i < numTasks - 1; i++){
    // Load sendbuffer with corpusbuffer
    // corpusbuffer to be renewed with MPI_Irecv
    // Now using sendbuffer for kNN calls since it is not altered during Isend
    double *temp;
    SWAP(sendbuffer, corpusbuffer);

    // No need to send anything on last iteration
    if(i != numTasks - 2){
      // Send received points to next process
      err = MPI_Isend(sendbuffer, size, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD, mpireq);
      if(err){
        MPI_Error_string(err, errbuffer, &errlen);
        printf("Error %d [%s]\n", err, errbuffer);
      }
      // Receive new points
      err = MPI_Irecv(corpusbuffer, size, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, mpireq+1);
      if(err){
        MPI_Error_string(err, errbuffer, &errlen);
        printf("Error %d [%s]\n", err, errbuffer);
      }
    }

    // Find kNN for new set of corpus points
    newknn = kNN(sendbuffer, X, n, n, d, k);
    // Offset indexes to match starting array
    for(j = 0; j < newknn.m * newknn.k; j++){
      newknn.nidx[j] += ((tid-2-i+numTasks)%numTasks) * n;
    }
    // and update selfknn for all corpus points checked so far
    mergekNN(&selfknn, &newknn);
    free(newknn.ndist);
    free(newknn.nidx);

    // Wait for sends/receives to finish before sending/recieving again
    MPI_Waitall(2, mpireq, mpistat);
  }

  free(corpusbuffer);
  free(sendbuffer);

  // Global minimum/maximum computation
  // double min = RAND_MAX, max = 0.0;
  // for(i = 0; i < selfknn.m * selfknn.k; i++){
  //   if(selfknn.ndist[i] < min)
  //     min = selfknn.ndist[i];
  //
  //   if(selfknn.ndist[i] > max)
  //     max = selfknn.ndist[i];
  // }
  //
  // printf("Local minimum: %f\n", min);
  // printf("Local maximum: %f\n", max);
  //
  // double maxAll, minAll;
  // MPI_Reduce(&max, &maxAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  // MPI_Reduce(&min, &minAll, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  //
  // if(tid == 0){
  //   printf("Global minimum: %f\n", minAll);
  //   printf("Global maximum: %f\n", maxAll);
  //   printf("\n");
  // }

  return selfknn;
}
