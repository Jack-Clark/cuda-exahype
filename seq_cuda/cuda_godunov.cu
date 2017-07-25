#include <stdio.h>
#include <math.h>
#include <string.h>
#include "kernels/KernelUtils.h"
#include "tarch/la/Vector.h"

#define DIMENSIONS 2


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__global__ void gen_eigenvalues_kernel(const double* const Q, const int normalNonZeroIndex, double* lambda, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const int ltid = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;
  int gtid;
  switch(idxSelector) {
    case 0: // This is the idx(j, k, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 1: // This is the idx(j, k+1, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + ((threadIdx.y+1) * numVariables);
      break;
    
    case 2: // This is the idx(j, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 3: // This is the idx(j+1, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return;
      }
      gtid = (threadIdx.x+1) * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    default:
      gtid = 0;
      assert(false);
      break;
  }

  // Application code goes here
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[gtid];
  const double p = (GAMMA-1) * (Q[gtid+4] - 0.5 * (Q[gtid+1] * Q[gtid+1] + Q[gtid+2] * Q[gtid+2]) * irho);   

  const double u_n = Q[gtid+normalNonZeroIndex + 1] * irho;
  const double c = sqrt(GAMMA * p * irho);

  lambda[ltid] = u_n - c;
  lambda[ltid+1] = u_n;
  lambda[ltid+2] = u_n;
  lambda[ltid+3] = u_n;
  lambda[ltid+4] = u_n + c;
}

__global__ void gen_flux_kernel_2D(const double* const Q, double** F, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const int ltid = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;
  int gtid = 0;
  switch(idxSelector) {
    case 0: // This is the idx(j, k, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 1: // This is the idx(j, k+1, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + ((threadIdx.y+1) * numVariables);
      break;
    
    case 2: // This is the idx(j, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return;
      }
      gtid = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 3: // This is the idx(j+1, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return;
      }
      gtid = (threadIdx.x+1) * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    default:
      gtid = 0;
      assert(false);
      break;
  }

  // Application code goes here
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[gtid];
  const double p = (GAMMA-1) * (Q[gtid+4] - 0.5 * (Q[gtid+1] * Q[gtid+1] + Q[gtid+2] * Q[gtid+2]) * irho);

  double* f = F[0];
  double* g = F[1];

  f[ltid] = Q[gtid+1];
  f[ltid+1] = irho * Q[gtid+1] * Q[gtid+1] + p;
  f[ltid+2] = irho * Q[gtid+1] * Q[gtid+2];
  f[ltid+3] = irho * Q[gtid+1] * Q[gtid+3];
  f[ltid+4] = irho * Q[gtid+1] * (Q[gtid+4] + p);

  g[ltid] = Q[gtid+2];
  g[ltid+1] = irho * Q[gtid+2] * Q[gtid+1];
  g[ltid+2] = irho * Q[gtid+2] * Q[gtid+2] + p;
  g[ltid+3] = irho * Q[gtid+2] * Q[gtid+3];
  g[ltid+4] = irho * Q[gtid+2] * (Q[gtid+4] + p);
}


__global__ void updateF_kernel(const double* const Q, double** FL2, double** FR2, double* F, const double* d_maxes, const int numVariables, 
                               const int patchBegin, const int basisSize, const int normalNonZero, const int xFlag) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  const int ltid = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;
  int qL_idx;
  int qR_idx;
  if(xFlag) {
      if(threadIdx.x < patchBegin) {
        return;
      }
      qL_idx = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      qR_idx = threadIdx.x * (numVariables * (basisSize+2)) + ((threadIdx.y+1) * numVariables);
  } else {
      if(threadIdx.y < patchBegin) {
        return;
      }
      qL_idx = threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      qR_idx = (threadIdx.x+1) * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
  }

  int tid = threadIdx.x + blockIdx.x * threadIdx.y;
  F[ltid]   = 0.5 * (FL2[normalNonZero][ltid]   + FR2[normalNonZero][ltid])   + 0.5 * d_maxes[tid] * (Q[qL_idx]   - Q[qR_idx]);
  F[ltid+1] = 0.5 * (FL2[normalNonZero][ltid+1] + FR2[normalNonZero][ltid+1]) + 0.5 * d_maxes[tid] * (Q[qL_idx+1] - Q[qR_idx+1]);
  F[ltid+2] = 0.5 * (FL2[normalNonZero][ltid+2] + FR2[normalNonZero][ltid+2]) + 0.5 * d_maxes[tid] * (Q[qL_idx+2] - Q[qR_idx+2]);
  F[ltid+3] = 0.5 * (FL2[normalNonZero][ltid+3] + FR2[normalNonZero][ltid+3]) + 0.5 * d_maxes[tid] * (Q[qL_idx+3] - Q[qR_idx+3]);
  F[ltid+4] = 0.5 * (FL2[normalNonZero][ltid+4] + FR2[normalNonZero][ltid+4]) + 0.5 * d_maxes[tid] * (Q[qL_idx+4] - Q[qR_idx+4]);
}

/* Finds the max value in sL and sR per cell and stores in global memory to be used later. */
__global__ void compute_maxes_kernel(double* sL, double* sR, double* maxes, const int numVariables) {
  
  const int ltid = threadIdx.x + blockIdx.x * threadIdx.y * numVariables;

  double max = -1.0;
  for(int i = 0; i < numVariables; i++) {
    const double abs_sL = fabs(sL[ltid+i]);
    const double abs_sR = fabs(sR[ltid+i]);
    const double tmp_max = fmax(abs_sL, abs_sR);
    max = fmax(max, tmp_max);
  }
  maxes[threadIdx.x + blockIdx.x * threadIdx.y] = max;

}

int c_index(int i, int j, int k, int basis, int numVariables) {
  return i * (basis * numVariables) + j * numVariables + k;
}

/**
 * Solves all the Riemann problems that do only require
 * internal data and add the result directly onto the
 * new solution. 
 * Finally add the source terms.
 */
extern "C"
double cudaSolutionUpdate(double* luh_new, const double* luh, 
                          int numberOfVariables, int basisSize,
                          double cflFactor, double dt_max_allowed,
                          const double cellSize[], int patchBegin,
                          int patchEnd, double dt) { 

  // Start CUDA
  int numCudaThreads = (patchEnd+1) * (patchEnd+1);
  int cudaMemReq = c_index(patchEnd+1, patchEnd, 5, basisSize+2, numberOfVariables);

  // Allocate and Transfer Memory
  double * d_Fn;
  double * d_sL;
  double * d_sR;
  double ** d_FL2;
  double ** d_FR2;
  double * d_q;
  gpuErrchk( cudaMalloc((void **) &d_Fn, numCudaThreads * numberOfVariables * sizeof(double)) );
  gpuErrchk( cudaMalloc((void **) &d_sL, numCudaThreads * numberOfVariables * sizeof(double)) );
  gpuErrchk( cudaMalloc((void **) &d_sR, numCudaThreads * numberOfVariables * sizeof(double)) );
  gpuErrchk( cudaMalloc((void **) &d_q, cudaMemReq * sizeof(double)) );

  gpuErrchk( cudaMemcpy(d_q, luh, cudaMemReq * sizeof(double), cudaMemcpyHostToDevice) );

  dim3 dimBlock(patchEnd+1, patchEnd+1);
  int normalNonZeroIndex = 0;

  // Compute x edges
  gen_eigenvalues_kernel<<<1, dimBlock>>>(d_q, normalNonZeroIndex, d_sL, numberOfVariables, patchBegin, basisSize, 0);
  gpuErrchk( cudaPeekAtLastError() );
  gen_eigenvalues_kernel<<<1, dimBlock>>>(d_q, normalNonZeroIndex, d_sR, numberOfVariables, patchBegin, basisSize, 1);
  gpuErrchk( cudaPeekAtLastError() );

  double * d_maxes;
  gpuErrchk( cudaMalloc((void **) &d_maxes, sizeof(double) * numberOfVariables * numCudaThreads) );
  cudaMemset(d_maxes, 0, sizeof(double) * numberOfVariables * numCudaThreads);

  compute_maxes_kernel<<<1, dimBlock>>>(d_sL, d_sR, d_maxes, numberOfVariables);
  gpuErrchk( cudaPeekAtLastError() );

  double maxes[numCudaThreads];
  memset(maxes, 0, sizeof(double) * numCudaThreads);

  gpuErrchk( cudaMemcpy(maxes, d_maxes, numCudaThreads * sizeof(double), cudaMemcpyDeviceToHost) );

  double s_max = -1.0;
  for (int j = 0; j < numCudaThreads; j++) {
    const double abs_temp = fabs(maxes[j]);
    s_max = fmax( abs_temp, s_max );
  }

  // TODO(guera): Improve. I'm quite sure this is not the correct/best
  // formula. TODO(Dominic): The division by DIMENSIONS might make sure that C_x+C_y < 1
  dt_max_allowed = fmin(
      dt_max_allowed, cflFactor / DIMENSIONS * cellSize[0] / s_max); // TODO(Dominic): Ignore this for a while
  double dt_max_1 = dt_max_allowed;

  int numDimensions = 2;
  /* This next section takes a bit of getting used to if you're not familiar with GPU programming.
     Basically, we create an array in host memory of device pointers that have been allocated device memory.
     This is necessary in order to work with pointer->pointer types in CUDA.
  */
  gpuErrchk( cudaMalloc((void **) &d_FL2, numDimensions * sizeof(double *)) );
  double* devicePointersStoredInHostMemoryL[numDimensions]; 
  for(int i = 0; i < numDimensions; i++) {
      gpuErrchk( cudaMalloc( (void**)&devicePointersStoredInHostMemoryL[i], numCudaThreads * numberOfVariables * sizeof(double)) );
  }
  gpuErrchk( cudaMemcpy(d_FL2, devicePointersStoredInHostMemoryL, sizeof(double*) * numDimensions, cudaMemcpyHostToDevice) );

  gpuErrchk( cudaMalloc((void **) &d_FR2, numDimensions * sizeof(double *)) );
  double* devicePointersStoredInHostMemoryR[numDimensions]; 
  for(int i = 0; i < numDimensions; i++) {
      gpuErrchk( cudaMalloc( (void**)&devicePointersStoredInHostMemoryR[i], numCudaThreads * numberOfVariables * sizeof(double)) );
  }
  gpuErrchk( cudaMemcpy(d_FR2, devicePointersStoredInHostMemoryR, sizeof(double*) * numDimensions, cudaMemcpyHostToDevice) );

  gen_flux_kernel_2D<<<1, dimBlock>>>(d_q, d_FL2, numberOfVariables, patchBegin, basisSize, 0);
  gpuErrchk( cudaPeekAtLastError() );
  gen_flux_kernel_2D<<<1, dimBlock>>>(d_q, d_FR2, numberOfVariables, patchBegin, basisSize, 1);
  gpuErrchk( cudaPeekAtLastError() );

  updateF_kernel<<<1, dimBlock>>>(d_q, d_FL2, d_FR2, d_Fn, d_maxes, numberOfVariables, patchBegin, basisSize, normalNonZeroIndex, 1);
  gpuErrchk( cudaPeekAtLastError() );

  double Fn[numCudaThreads * numberOfVariables];
  gpuErrchk( cudaMemcpy(Fn, d_Fn, numCudaThreads * numberOfVariables * sizeof(double), cudaMemcpyDeviceToHost) );
  
  for (int j = patchBegin; j < patchEnd+1; j++) {
    for (int k = patchBegin-1; k < patchEnd+1; k++) {
      for (int l=0; l<numberOfVariables; ++l) {
        luh_new[c_index(j, k, l, basisSize+2, numberOfVariables)]   -= dt / cellSize[0] * Fn[l];  
        luh_new[c_index(j, k+1, l, basisSize+2, numberOfVariables)] += dt / cellSize[0] * Fn[l];
      }
    }
  }

  // Compute y faces
  normalNonZeroIndex = 1;
  gen_eigenvalues_kernel<<<1, dimBlock>>>(d_q, normalNonZeroIndex, d_sL, numberOfVariables, patchBegin, basisSize, 2);
  gpuErrchk( cudaPeekAtLastError() );
  gen_eigenvalues_kernel<<<1, dimBlock>>>(d_q, normalNonZeroIndex, d_sR, numberOfVariables, patchBegin, basisSize, 3);
  gpuErrchk( cudaPeekAtLastError() );

  compute_maxes_kernel<<<1, dimBlock>>>(d_sL, d_sR, d_maxes, numberOfVariables);
  gpuErrchk( cudaPeekAtLastError() );

  gpuErrchk( cudaMemcpy(maxes, d_maxes, numCudaThreads * sizeof(double), cudaMemcpyDeviceToHost) );

  s_max = -1.0;
  for (int j = 0; j < (patchEnd+1)*(patchEnd+1); j++) {
    const double abs_temp = fabs(maxes[j]);
    //printf("max[%d] = %f\n", j, abs_temp);
    s_max = fmax( abs_temp, s_max );
  }

  // TODO(guera): Improve. I'm quite sure this is not the correct/best
  // formula. TODO(Dominic): The division by DIMENSIONS might make sure that C_x+C_y < 1
  dt_max_allowed = fmin(
      dt_max_allowed, cflFactor / DIMENSIONS * cellSize[1] / s_max); // TODO(Dominic): Ignore this for a while
  double dt_max_2 = dt_max_allowed;

  gen_flux_kernel_2D<<<1, dimBlock>>>(d_q, d_FL2, numberOfVariables, patchBegin, basisSize, 2);
  gpuErrchk( cudaPeekAtLastError() );
  gen_flux_kernel_2D<<<1, dimBlock>>>(d_q, d_FR2, numberOfVariables, patchBegin, basisSize, 3);
  gpuErrchk( cudaPeekAtLastError() );

  updateF_kernel<<<1, dimBlock>>>(d_q, d_FL2, d_FR2, d_Fn, d_maxes, numberOfVariables, patchBegin, basisSize, normalNonZeroIndex, 0);
  gpuErrchk( cudaPeekAtLastError() );

  gpuErrchk( cudaMemcpy(Fn, d_Fn, numCudaThreads * numberOfVariables * sizeof(double), cudaMemcpyDeviceToHost) );
  
  for (int j = patchBegin-1; j < patchEnd+1; j++) {
    for (int k = patchBegin; k < patchEnd+1; k++) {
      for (int l=0; l<numberOfVariables; ++l) {
        luh_new[c_index(j, k, l, basisSize+2, numberOfVariables)]   -= dt / cellSize[1] * Fn[l];  
        luh_new[c_index(j+1, k, l, basisSize+2, numberOfVariables)] += dt / cellSize[1] * Fn[l];
      }
    }
  }

  cudaFree(d_Fn);
  cudaFree(d_sL);
  cudaFree(d_sR);
  cudaFree(d_maxes);
  cudaFree(d_q);
  for(int i = 0; i < DIMENSIONS; i++) {
      cudaFree(devicePointersStoredInHostMemoryL[i]);
      cudaFree(devicePointersStoredInHostMemoryR[i]);
  }
  cudaFree(d_FR2);
  cudaFree(d_FL2);
    
  return dt_max_allowed;
}
