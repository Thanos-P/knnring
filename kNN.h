/**
 * @file   kNN.h
 * @author athanasps <athanasps@ece.auth.gr>
 *         Thanos Paraskevas
 *
 */

#include <stdio.h>
#include <stdlib.h>


// Definition of the kNN result struct
typedef struct knnresult{
int * nidx;           //!< Indices (0-based) of nearest neighbors   [m-by-k]
double * ndist;       //!< Distance of nearest neighbors            [m-by-k]
int m;                //!< Number of query points                   [scalar]
int k;                //!< Number of nearest neighbors              [scalar]
} knnresult;


//! Compute k nearest neighbors of each point in X [n-by-d]
/*!

\param X Corpus data points             [n-by-d]
\param Y Query data points              [m-by-d]
\param n Number of data points          [scalar]
\param d Number of dimensions           [scalar]
\param k Number of neighbors            [scalar]

\return The kNN result
*/
knnresult kNN(double * X, double * Y, int n, int d, int k);
