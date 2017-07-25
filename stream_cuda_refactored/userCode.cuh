#ifndef __CUDA_USER_CODE_CUH__
#define __CUDA_USER_CODE_CUH__

__global__ void eigenvalues_kernel(const double* const Q, const int normalNonZeroIndex, double* lambda, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector);

__global__ void flux_kernel(const double* const Q, double** F, const int numVariables, const int patchBegin, const int basisSize, const int idxSelector);

#endif
