/**
 * @file   kNN_sequential.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cblas.h"
#include "cblas_f77.h"
#include "kNN.h"

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
\param d Number of dimensions           [scalar]
\param k Number of neighbors            [scalar]

\return The kNN result
*/
knnresult kNN(double * X, double * Y, int n, int m, int d, int k) {
  int i, j;
  double *dist = (double *)malloc(m*n*sizeof(double));
  int *idx = (int *)malloc(m*n*sizeof(int));

  // Find distances

  int inc = 1;
  for(i = 0; i < n; i++){
    // MATLAB: sum(X.^2,2)
    dist[i*m] = pow(cblas_dnrm2(d, (X+i*d), inc), 2);
    for(j = 0; j < m; j++){
      // Value to be expanded
      double row_expand;
      if(j == 0)
        row_expand = dist[i*m];
      // Expand value for all m and add MATLAB: sum(Y.^2,2).'
      dist[i*m + j] = row_expand + pow(cblas_dnrm2(d, (Y+j*d), inc), 2);
    }
  }

  enum CBLAS_ORDER order = CblasRowMajor;
  enum CBLAS_TRANSPOSE transX = CblasNoTrans;
  enum CBLAS_TRANSPOSE transY = CblasTrans;
  double alpha = -2.0, beta = 1.0;
  int ldX = d, ldY = d, lddist = m;

  // MATLAB: - 2 * X*Y.' + sum(X.^2,2) + sum(Y.^2,2).'
  cblas_dgemm(order, transX, transY, n, m, d, alpha, X, ldX, Y, ldY, beta, dist, lddist);
  // MATLAB: sqrt(ans)
  for(i = 0; i < n*m; i++)
    dist[i] = sqrt(dist[i]);

  double *distT = (double *)malloc(m*n*sizeof(double));

  for(i = 0; i < n; i++){
    for(j = 0; j < m; j++){
      distT[j*n + i] = dist[i*m + j];
      idx[j*n + i] = i;
    }
  }

  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      printf("(%d) %f ", idx[i*n + j], distT[i*n + j]);
    }
    printf("\n");
  }

  // Initialize result
  knnresult result;
  result.m = 1;
  result.k = k;
  result.nidx = (int *)malloc(m*k*sizeof(int));
  result.ndist = (double *)malloc(m*k*sizeof(double));

  // Find k nearest neighbors using quickselect for each query point
  for(i = 0; i < m; i++){
    printf("\nNearest neighbors of query point %d:\n", i);
    for(j = 0; j < k; j++){
      quickselect((distT+i*n), (idx+i*n), 0, n-1, j);
      result.nidx[i*k + j] = idx[i*n + j];
      result.ndist[i*k + j] = distT[i*n + j];
      printf("%d) index = %d, distance = %f\n", j, idx[i*n + j], distT[i*n + j]);
    }
  }

  free(distT);
  free(dist);
  free(idx);

  return result;
}
