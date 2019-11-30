/**
 * @file   main_mpi.c
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <mpi.h>
#include <sys/time.h>
#include "knnring.h"


double mpikNN(int n, int d, int k) {

  int numTasks, tid;
  double *corpus = NULL;
  double *corpusAll = NULL;

  MPI_Comm_rank(MPI_COMM_WORLD, &tid);
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);

  if(tid == 0){                // ..... MASTER

    corpusAll = (double *)malloc(numTasks*n*d * sizeof(double));

    for (int ip = 0; ip < numTasks; ip++){
      // "read" new chunk
      corpus = (double *)malloc(n*d * sizeof(double));

      for (int i=0;i<n*d;i++){
        corpusAll[i+ip*n*d] = (double) 10 * rand() / RAND_MAX;
        corpus[i]= corpusAll[i+ip*n*d];
      }

      if(ip == numTasks-1)            // last chunk is mine
        break;

      int dst = ip+1;
      int tag = 1;
      MPI_Send(corpus, n*d, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);

      free(corpus);
    }

  } else {                      // ..... other processes

    int rcv = 0;
    int tag = 1;
    MPI_Status Stat;
    corpus = (double *)malloc(n*d * sizeof(double));
    MPI_Recv(corpus, n*d, MPI_DOUBLE, rcv, tag, MPI_COMM_WORLD, &Stat);
  }

  struct timeval startwtime, endwtime;
  double totaltime;

  // run the distributed kNN code and time it
  gettimeofday(&startwtime, NULL);

  knnresult knnres = distrAllkNN(corpus, n, d, k);

  gettimeofday(&endwtime, NULL);

  totaltime = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
            + endwtime.tv_sec - startwtime.tv_sec);

  // ~~~~~~~~~~~~~~~~~~~~ gather back kNN results

  if (tid == 0) {                // ..... MASTER

    knnresult knnresall;
    knnresall.nidx  = (int *)   malloc( n*numTasks*k*sizeof(int)    );
    knnresall.ndist = (double *)malloc( n*numTasks*k*sizeof(double) );
    knnresall.m = n*numTasks;
    knnresall.k = k;

    for (int ip = 0; ip < numTasks-1; ip++){
      int rcv = ip+1;
      int tag = 1;
      MPI_Status Stat;

      MPI_Recv(&knnresall.nidx[ip*n*k], n*k, MPI_INT, rcv, tag,
                MPI_COMM_WORLD, &Stat);

      MPI_Recv(&knnresall.ndist[ip*n*k], n*k, MPI_DOUBLE, rcv, tag,
                MPI_COMM_WORLD, &Stat);
    }

    // move my result to final struct
    for (int i = 0; i < n*k; i++){
      knnresall.nidx[(numTasks-1)*n*k+i] = knnres.nidx[i];
      knnresall.ndist[(numTasks-1)*n*k+i] = knnres.ndist[i];
    }

    // Test print for overall results
    // int i, j;
    // for(i = 0; i < knnresall.m; i++){
    //   printf("\nNearest neighbors of query point %d:\n", i);
    //   for(j = 0; j < knnresall.k; j++){
    //     printf("%d) index = %d, distance = %f\n",
    //       j, knnresall.nidx[i*k + j], knnresall.ndist[i*k + j]);
    //   }
    // }

    free(knnresall.nidx);
    free(knnresall.ndist);

  } else {                      // ..... other processes
      int dst = 0;
      int tag = 1;

      // send to correct process
      MPI_Send(knnres.nidx, n*k, MPI_INT, dst, tag, MPI_COMM_WORLD);
      MPI_Send(knnres.ndist, n*k, MPI_DOUBLE, dst, tag, MPI_COMM_WORLD);
  }

  free(corpus);
  free(corpusAll);

  return totaltime;
}

int main(int argc, char *argv[]){

  MPI_Init(&argc, &argv);

  int n=1423;       // corpus elements per process
  int d=37;         // dimensions
  int k=13;         // nearest neighbors

  int numTasks, tid;
  double totaltime;

  MPI_Comm_rank(MPI_COMM_WORLD, &tid);
  MPI_Comm_size(MPI_COMM_WORLD, &numTasks);

  totaltime = mpikNN(n, d, k);

  printf("Process: %d / %d - Time: %f\n", tid, numTasks - 1, totaltime);

  MPI_Barrier(MPI_COMM_WORLD);

  // Find mean time
  if(tid == 0){
    MPI_Status mpistat;
    double otherstime;
    for(int i = 1; i < numTasks; i++){
      MPI_Recv(&otherstime, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &mpistat);
      totaltime += otherstime;
    }
    printf("Mean time: %f\n", totaltime/numTasks);
  }else{
    MPI_Send(&totaltime, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
  }

  MPI_Finalize();

  return 0;
}
