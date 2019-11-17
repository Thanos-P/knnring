/**
 * @file   kNN_sequential.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "knnring.h"

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

//! Compute k nearest neighbors of each point in X [n-by-d]
/*!

\param X Corpus data points             [n-by-d]
\param Y Query data points              [m-by-d]
\param n Number of data points          [scalar]
\param m Number of query points         [scalar]
\param d Number of dimensions           [scalar]
\param k Number of neighbors            [scalar]

\return The kNN result
*/
knnresult kNN(double * X, double * Y, int n, int m, int d, int k) {
  int i, j;
  // Initialize distances and indexes arrays
  double *dist = (double *)malloc(m*n*sizeof(double));
  int *idx = (int *)malloc(m*n*sizeof(int));

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
  // MATLAB: sqrt(ans)
  for(i = 0; i < n*m; i++)
    dist[i] = sqrt(dist[i]);

  // Test print
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      printf("(%d) %f ", idx[i*n + j], dist[i*n + j]);
    }
    printf("\n\n");
  }

  // Initialize result
  knnresult result;
  result.m = m;
  result.k = k;
  result.nidx = (int *)malloc(m*k*sizeof(int));
  result.ndist = (double *)malloc(m*k*sizeof(double));

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
