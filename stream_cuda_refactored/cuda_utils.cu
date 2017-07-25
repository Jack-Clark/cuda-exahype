#include <stdio.h>

void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__device__ int idx(int threadIdxX, int threadIdxY, int blockIdxX, int numVariables, int basisSize, int patchBegin, int idxSelector) {
  switch(idxSelector) {
    case 0: // This is the idx(j, k, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return -1;
      }
      return threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 1: // This is the idx(j, k+1, 0) index for the x face
      if(threadIdx.x < patchBegin) {
        return -1;
      }
      return threadIdx.x * (numVariables * (basisSize+2)) + ((threadIdx.y+1) * numVariables);
      break;
    
    case 2: // This is the idx(j, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return -1;
      }
      return threadIdx.x * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    case 3: // This is the idx(j+1, k, 0) index for the y edge
      if(threadIdx.y < patchBegin) {
        return -1;
      }
      return (threadIdx.x+1) * (numVariables * (basisSize+2)) + (threadIdx.y * numVariables);
      break;

    default:
      return -1;
      break;
  }
}