/**
 * @file   kNN_sequential.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "kNN.h"

//! Swaps the values of element a and b
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

//! Return the index of k-th smallest element of array
//  within left..right inclusive
void quickselect(double *array, int *idx, int left, int right, int k){
  if(left == right){      // If the list contains only one element,
    return;     // return that element
  }
  // select a random pivotIndex between left and right
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

knnresult kNN(double * X, double * Y, int n, int d, int k) {
  int i, j;
  double *dist = (double *)malloc(n*sizeof(double));
  int *idx = (int *)malloc(n*sizeof(int));
  printf("Distances:\n");
  for(i = 0; i < n; i++){
    idx[i] = i;
    dist[i] = 0;
    for(j = 0; j < d; j++){
      dist[i] += pow(X[i*d + j] - Y[j], 2);
    }
    dist[i] = sqrt(dist[i]);
    printf("%d) %f\n", i, dist[i]);
  }

  knnresult result;
  result.m = 1;
  result.k = k;
  result.nidx = (int *)malloc(k*sizeof(int));
  result.ndist = (double *)malloc(k*sizeof(double));
  printf("\nNearest neighbors:\n");
  for(i = 0; i < k; i++){
    quickselect(dist, idx, 0, n-1, i);
    result.nidx[i] = idx[i];
    result.ndist[i] = dist[i];
    printf("%d) index = %d, distance = %f\n", i, idx[i], dist[i]);
  }

  free(dist);
  free(idx);

  return result;
}
