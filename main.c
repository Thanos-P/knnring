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
#include "knnring.h"

double * createRandomPoints(int n, int d){

  int i, j;
  // Initialize points array
  double *X = (double *)malloc(n*d*sizeof(double));
  printf("POINTS\n");
  for(i = 0; i < n; i++){
    printf("%d) ( ", i);
    for(j = 0; j < d; j++){
      *(X + i*d + j) = (double) 10 * rand() / RAND_MAX;
      printf("%f, ", *(X + i*d + j));
    }
    printf(")\n");
  }
  printf("\n");
  return X;
}

int main(){

  int n = 16;
  int m = 16;
  int d = 2;
  int k = 4;

  srand(time(NULL));
  double *X = createRandomPoints(n, d);
  double *Y = createRandomPoints(m, d);


  knnresult result = kNN(X, X, n, n, d, k);

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
