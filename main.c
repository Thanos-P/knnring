/**
 * @file   main.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "knnring.h"

double * createRandomPoints(int n, int d){

  int i, j;
  // Initialize points array
  double *X = (double *)malloc(n*d*sizeof(double));
  for(i = 0; i < n; i++){
    for(j = 0; j < d; j++){
      *(X + i*d + j) = (double) 10 * rand() / RAND_MAX;
    }
  }
  return X;
}

int main(){

  int n = 1423;
  int m = 1;
  int d = 37;
  int k = 13;

  srand(time(NULL));
  double *X = createRandomPoints(n, d);
  double *Y = createRandomPoints(m, d);

  struct timeval startwtime, endwtime;
  double totaltime;

  gettimeofday(&startwtime, NULL);

  knnresult result = kNN(X, X, n, n, d, k);

  gettimeofday(&endwtime, NULL);

  totaltime = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
            + endwtime.tv_sec - startwtime.tv_sec);

  printf("Total time: %f\n", totaltime);

  // Test print for overall results
  // int i, j;
  // for(i = 0; i < result.m; i++){
  //   printf("\nNearest neighbors of query point %d:\n", i);
  //   for(j = 0; j < result.k; j++){
  //     printf("%d) index = %d, distance = %f\n",
  //      j, result.nidx[i*k + j], result.ndist[i*k + j]);
  //   }
  // }

  free(result.nidx);
  free(result.ndist);
  free(X);
  free(Y);

  return 0;
}
