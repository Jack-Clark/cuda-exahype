#include <stdio.h>
#include <math.h>
#include "cuda_utils.cuh"

/* 
 * The user's eigenvalues CUDA implementation should go here. All that is required to change from a CPU version to a CUDA version
 * is for the user to add the value of q_arr_idx to their Q array indices and to add other_idx to all other array indices. For example:
 * Q[1] becomes Q[q_arr_idx + 1] and lambda[1] becomes lambda[other_idx + 1].
 */
__global__ void eigenvalues_kernel(const double* const Q, const int normalNonZeroIndex, double* lambda, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const int other_idx = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;
  const int q_arr_idx = idx(threadIdx.x, threadIdx.y, blockIdx.x, numVariables, basisSize, patchBegin, idxSelector);
  if(q_arr_idx == -1) 
    return;

  // Application code goes here
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[q_arr_idx];
  const double p = (GAMMA-1) * (Q[q_arr_idx+4] - 0.5 * (Q[q_arr_idx+1] * Q[q_arr_idx+1] + Q[q_arr_idx+2] * Q[q_arr_idx+2]) * irho);   

  const double u_n = Q[q_arr_idx+normalNonZeroIndex + 1] * irho;
  const double c = sqrt(GAMMA * p * irho);

  lambda[other_idx] = u_n - c;
  lambda[other_idx+1] = u_n;
  lambda[other_idx+2] = u_n;
  lambda[other_idx+3] = u_n;
  lambda[other_idx+4] = u_n + c;
}

/* 
 * The user's flux CUDA implementation should go here. All that is required to change from a CPU version to a CUDA version
 * is for the user to add the value of q_arr_idx to their Q array indices and to add other_idx to all other array indices. For example:
 * Q[1] becomes Q[q_arr_idx + 1] and f[1] becomes f[other_idx + 1].
 */
__global__ void flux_kernel(const double* const Q, double** F, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const int other_idx = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;
  const int q_arr_idx = idx(threadIdx.x, threadIdx.y, blockIdx.x, numVariables, basisSize, patchBegin, idxSelector);
  if(q_arr_idx == -1) 
    return;

  const double GAMMA = 1.4;

  const double irho = 1.0/Q[q_arr_idx];
  const double p = (GAMMA-1) * (Q[q_arr_idx+4] - 0.5 * (Q[q_arr_idx+1] * Q[q_arr_idx+1] + Q[q_arr_idx+2] * Q[q_arr_idx+2]) * irho);

  double* f = F[0];
  double* g = F[1];

  f[other_idx] = Q[q_arr_idx+1]; // should be numVariables * tid
  f[other_idx+1] = irho * Q[q_arr_idx+1] * Q[q_arr_idx+1] + p;
  f[other_idx+2] = irho * Q[q_arr_idx+1] * Q[q_arr_idx+2];
  f[other_idx+3] = irho * Q[q_arr_idx+1] * Q[q_arr_idx+3];
  f[other_idx+4] = irho * Q[q_arr_idx+1] * (Q[q_arr_idx+4] + p);

  g[other_idx] = Q[q_arr_idx+2]; // Should be numVariables * tid + numVariables
  g[other_idx+1] = irho * Q[q_arr_idx+2] * Q[q_arr_idx+1];
  g[other_idx+2] = irho * Q[q_arr_idx+2] * Q[q_arr_idx+2] + p;
  g[other_idx+3] = irho * Q[q_arr_idx+2] * Q[q_arr_idx+3];
  g[other_idx+4] = irho * Q[q_arr_idx+2] * (Q[q_arr_idx+4] + p);
}