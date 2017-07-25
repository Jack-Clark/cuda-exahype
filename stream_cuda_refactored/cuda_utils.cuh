#ifndef __CUDA_UTILS_CUH__
#define __CUDA_UTILS_CUH__

#define DIMENSIONS 2

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true);

__device__ int idx(int threadIdxX, int threadIdxY, int blockIdxX, int numVariables, int basisSize, int patchBegin, int idxSelector);

#endif