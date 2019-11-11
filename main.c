#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "kNN.h"

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
  int n = 8, m = 4, d = 2, k = 2;
  srand(time(NULL));
  double *X = createRandomPoints(n, d);
  double *Y = createRandomPoints(m, d);
  knnresult result = kNN(X, Y, n, m, d, k);

  free(result.nidx);
  free(result.ndist);
  free(X);
  free(Y);

  return 0;
}
