/**
 * This file is part of the ExaHyPE project.
 * Copyright (c) 2016  http://exahype.eu
 * All rights reserved.
 *
 * The project has received funding from the European Union's Horizon
 * 2020 research and innovation programme under grant agreement
 * No 671698. For copyrights and licensing, please consult the webpage.
 *
 * Released under the BSD 3 Open Source License.
 * For the full license text, see LICENSE.txt
 **/

#include "../Limiter.h"

#include "tarch/la/ScalarOperations.h"

namespace kernels {
namespace limiter {
namespace generic {
namespace c {
  
//Fortran (Limiter.f90): GetSubcellData
// Allocate lim memory via
// double* lim = new double[basisSizeLim*basisSizeLim*basisSizeLim*numberOfVariables]; //Fortran ref: lim(nVar,nSubLimV(1),nSubLimV(2),nSubLimV(3))
void projectOnFVLimiterSpace(const double* const luh, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const lim) {
  const int basisSizeLim = getBasisSizeLim(basisSize);
  
#if DIMENSIONS == 3
  const int basisSize3D       = basisSize;
  const int basisSizeLim3D    = basisSizeLim;
  const int ghostLayerWidth3D = ghostLayerWidth;
#else
  constexpr int basisSize3D       = 1;
  constexpr int basisSizeLim3D    = 1;
  constexpr int ghostLayerWidth3D = 0;
#endif

  idx4 idxLuh(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx4 idxLim(basisSizeLim3D+2*ghostLayerWidth3D, basisSizeLim+2*ghostLayerWidth, basisSizeLim+2*ghostLayerWidth, numberOfVariables);
  idx2 idxConv(basisSize, basisSizeLim); 
  
  int x,y,z,v,ix,iy,iz; 
  
  //tensor operation
  for(z=ghostLayerWidth3D; z<basisSizeLim3D+ghostLayerWidth3D; z++) { // We can skip x,y,z>=basisSizeLim+ghostLayerWidth
    for(y=ghostLayerWidth; y<basisSizeLim+ghostLayerWidth; y++) {
      for(x=ghostLayerWidth; x<basisSizeLim+ghostLayerWidth; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lim[idxLim(z,y,x,v)] = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize; ix++) {
                lim[idxLim(z,y,x,v)] += luh[idxLuh(iz,iy,ix,v)] 
                                        #if DIMENSIONS == 3
                                        * uh2lim[basisSize-1][idxConv(iz,z-ghostLayerWidth3D)]
                                        #endif
                                        * uh2lim[basisSize-1][idxConv(iy,y-ghostLayerWidth)]
                                        * uh2lim[basisSize-1][idxConv(ix,x-ghostLayerWidth)];
              }
            }
          }
        }
      }
    }
  }
  
}

//Fortran (Limiter.f90): PutSubcellData
void projectOnDGSpace(const double* const lim, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const luh) {
  const int basisSizeLim = getBasisSizeLim(basisSize);
  
#if DIMENSIONS == 3
  const int basisSize3D       = basisSize;
  const int basisSizeLim3D    = basisSizeLim;
  const int ghostLayerWidth3D = ghostLayerWidth;
#else
  constexpr int basisSize3D       = 1;
  constexpr int basisSizeLim3D    = 1;
  constexpr int ghostLayerWidth3D = 0;
#endif
  
  idx4 idxLuh(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx4 idxLim(basisSizeLim3D+2*ghostLayerWidth3D, basisSizeLim+2*ghostLayerWidth, basisSizeLim+2*ghostLayerWidth, numberOfVariables);
  idx2 idxConv(basisSizeLim, basisSize);
  
  int x,y,z,v,ix,iy,iz;

  //tensor operation
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          luh[idxLuh(z,y,x,v)] = 0;
          for(iz=ghostLayerWidth3D; iz<basisSizeLim3D+ghostLayerWidth3D; iz++) { // We can skip ix,iy,iz>=basisSizeLim+ghostLayerWidth
            for(iy=ghostLayerWidth; iy<basisSizeLim+ghostLayerWidth; iy++) {
              for(ix=ghostLayerWidth; ix<basisSizeLim+ghostLayerWidth; ix++) {
                luh[idxLuh(z,y,x,v)] += lim[idxLim(iz,iy,ix,v)] 
                                        #if DIMENSIONS == 3
                                        * lim2uh[basisSize-1][idxConv(iz-ghostLayerWidth3D,z)]
                                        #endif
                                        * lim2uh[basisSize-1][idxConv(iy-ghostLayerWidth,y)]
                                        * lim2uh[basisSize-1][idxConv(ix-ghostLayerWidth,x)];
              }
            }
          }
        }
      }
    }
  }
  
}

bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, const int numberOfVariables, const int basisSize,
                    const double DMPMaximumRelaxationParameter,const double DMPDifferenceScaling,
                    double* boundaryMinPerVariables, double* boundaryMaxPerVariables) {
  const int order = basisSize-1;
#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;
#endif

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  bool discreteMaximumPrincipleIsSatisfied;
  int x, y, z, iVar ,ix, iy, iz;
  double newLocalMin, newLocalMax, lobValue, boundaryMin, boundaryMax, ldiff;

  discreteMaximumPrincipleIsSatisfied = true;
  for(iVar = 0; iVar < numberOfVariables; iVar++) {
    newLocalMin = luh[iVar];
    newLocalMax = newLocalMin;

    //get min/max from luh
    for(z=0; z<basisSize3D; z++) {
      for(y=0; y<basisSize; y++) {
        for(x=0; x<basisSize; x++) {
          newLocalMin = std::min ( newLocalMin, luh[idx(z, y, x, iVar)] );
          newLocalMax = std::max ( newLocalMax, luh[idx(z, y, x, iVar)] );
        }
      }
    }

    //get min/max from lob: project and compare
    for(z=0; z<basisSize3D; z++) {
      for(y=0; y<basisSize; y++) {
        for(x=0; x<basisSize; x++) {
          lobValue = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                  lobValue += luh[idx(iz,iy,ix,iVar)]
                              #if DIMENSIONS == 3
                              * uh2lob[order][idxConv(iz,z)]
                              #endif
                              * uh2lob[order][idxConv(iy,y)] * uh2lob[order][idxConv(ix,x)];
              }
            }
          }
          newLocalMin = std::min( newLocalMin, lobValue );
          newLocalMax = std::max( newLocalMax, lobValue );
        }
      }
    }

    //evaluate troubled status for the given iVar
    boundaryMin = boundaryMinPerVariables[iVar];
    for (int x=1; x<DIMENSIONS_TIMES_TWO; x++) {
      boundaryMin = std::min( boundaryMin, boundaryMinPerVariables[x*numberOfVariables+iVar] );
    }
    boundaryMax = boundaryMaxPerVariables[iVar];
    for (int x=1; x<DIMENSIONS_TIMES_TWO; x++) {
      boundaryMax = std::max( boundaryMax, boundaryMaxPerVariables[x*numberOfVariables+iVar] );
    }
    ldiff = (boundaryMax - boundaryMin) * DMPDifferenceScaling;
    assertion3(tarch::la::greaterEquals(ldiff,0.0),ldiff,boundaryMin,boundaryMax);

    ldiff = std::max( ldiff, DMPMaximumRelaxationParameter );

    if((newLocalMin < (boundaryMin - ldiff)) ||
       (newLocalMax > (boundaryMax + ldiff))) {
      discreteMaximumPrincipleIsSatisfied = false;
    }

    // We have the new min and max directly available now.
    boundaryMinPerVariables[iVar] = newLocalMin;
    boundaryMaxPerVariables[iVar] = newLocalMax;
  }

  return discreteMaximumPrincipleIsSatisfied;
}


bool discreteMaximumPrinciple(const double* const luh, const int numberOfVariables, const int basisSize,
                    const double DMPMaximumRelaxationParameter,const double DMPDifferenceScaling,
                    const double* const boundaryMinPerVariables, const double* const boundaryMaxPerVariables) {
  const int order = basisSize-1;
#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;
#endif

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  int x, y, z, iVar ,ix, iy, iz;
  double min, max, lobValue, boundaryMin, boundaryMax, ldiff;

  for(iVar = 0; iVar < numberOfVariables; iVar++) {
    min = luh[iVar];
    max = min;

    //get min/max from luh
    for(z=0; z<basisSize3D; z++) {
      for(y=0; y<basisSize; y++) {
        for(x=0; x<basisSize; x++) {
          min = std::min ( min, luh[idx(z, y, x, iVar)] );
          max = std::max ( max, luh[idx(z, y, x, iVar)] );
        }
      }
    }

    //get min/max from lob: project and compare
    for(z=0; z<basisSize3D; z++) {
      for(y=0; y<basisSize; y++) {
        for(x=0; x<basisSize; x++) {
          lobValue = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                  lobValue += luh[idx(iz,iy,ix,iVar)]
                              #if DIMENSIONS == 3
                              * uh2lob[order][idxConv(iz,z)]
                              #endif
                              * uh2lob[order][idxConv(iy,y)] * uh2lob[order][idxConv(ix,x)];
              }
            }
          }
          min = std::min( min, lobValue );
          max = std::max( max, lobValue );
        }
      }
    }

    //evaluate troubled status for the given iVar
    boundaryMin = boundaryMinPerVariables[iVar];
    for (int x=1; x<DIMENSIONS_TIMES_TWO; x++) {
      boundaryMin = std::min( boundaryMin, boundaryMinPerVariables[x*numberOfVariables+iVar] );
    }
    boundaryMax = boundaryMaxPerVariables[iVar];
    for (int x=1; x<DIMENSIONS_TIMES_TWO; x++) {
      boundaryMax = std::max( boundaryMax, boundaryMaxPerVariables[x*numberOfVariables+iVar] );
    }
    ldiff = (boundaryMax - boundaryMin) * DMPDifferenceScaling;
    assertion3(tarch::la::greaterEquals(ldiff,0.0),ldiff,boundaryMin,boundaryMax);

    ldiff = std::max( ldiff, DMPMaximumRelaxationParameter );

    if((min < (boundaryMin - ldiff)) ||
       (max > (boundaryMax + ldiff))) {
      return false;
    }
  }

  return true;
}

/**
 * localMinPerVariables, localMaxPerVariables are double[numberOfVariables]
 */
void findCellLocalMinAndMax(const double* const luh, const int numberOfVariables, const int basisSize, double* const localMinPerVariables, double* const localMaxPerVariables) {
  int index, ii, iVar, iiEnd;
  
  // initialize and process luh
  index = 0;
  iiEnd =  basisSize*basisSize;
#if DIMENSIONS == 3
  iiEnd *= basisSize;
#endif
  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMinPerVariables[iVar] = luh[index];
    localMaxPerVariables[iVar] = luh[index];   
    index++;
  }
  for(ii = 1; ii < iiEnd; ii++) {
    for(iVar = 0; iVar < numberOfVariables; iVar++) {
      localMinPerVariables[iVar] = std::min ( localMinPerVariables[iVar], luh[index] );
      localMaxPerVariables[iVar] = std::max ( localMaxPerVariables[iVar], luh[index] );
      index++;
    }
  }
 
  //process lob
  compareWithADERDGSolutionAtGaussLobattoNodes(luh, numberOfVariables, basisSize, localMinPerVariables, localMaxPerVariables);
}

/**
 * localMinPerVariables, localMaxPerVariables are double[numberOfVariables]
 *
 * basisSize is the basis size of the ADER-DG scheme.
 */
void findCellLocalLimiterMinAndMax(const double* const lim, const int numberOfVariables, const int basisSize, const int ghostLayerWidth,
                                   double* const localMinPerVariables, double* const localMaxPerVariables) {
  const int basisSizeLim = getBasisSizeLim(basisSize);

#if DIMENSIONS == 3
  const int basisSizeLim3D    = basisSizeLim;
  const int ghostLayerWidth3D = ghostLayerWidth;
#else
  constexpr int basisSizeLim3D    = 1;
  constexpr int ghostLayerWidth3D = 0;
#endif
  idx4 idxLim(basisSizeLim3D+2*ghostLayerWidth3D, basisSizeLim+2*ghostLayerWidth, basisSizeLim+2*ghostLayerWidth, numberOfVariables);

  for(int iVar = 0; iVar < numberOfVariables; iVar++) {
    localMinPerVariables[iVar] = std::numeric_limits<double>::max();
    localMaxPerVariables[iVar] = -std::numeric_limits<double>::max();
  }

  for (int k=ghostLayerWidth3D; k<basisSizeLim3D+ghostLayerWidth3D; ++k) // skip the last element
  for (int j=ghostLayerWidth; j<basisSizeLim+ghostLayerWidth; ++j)
  for (int i=ghostLayerWidth; i<basisSizeLim+ghostLayerWidth; ++i) {
    for(int iVar = 0; iVar < numberOfVariables; iVar++) {
      localMinPerVariables[iVar] = std::min ( localMinPerVariables[iVar], lim[idxLim(k,j,i,iVar)] );
      localMaxPerVariables[iVar] = std::max ( localMaxPerVariables[iVar], lim[idxLim(k,j,i,iVar)] );
    }
  }
}

//*************************
//*** Private functions ***
//*************************

//Used only for test purpose to check the projection since compareWithADERDGSolutionAtGaussLobattoNodes and isTroubled uses the same loop
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize) {

#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;  
#endif

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double* lob = new double[basisSize3D*basisSize*basisSize*numberOfVariables]; //Fortran ref: lob(nVar,nDOF(1),nDOF(2),nDOF(3))
  
  int x,y,z,v,ix,iy,iz;
  
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lob[idx(z,y,x,v)] = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                lob[idx(z,y,x,v)] += luh[idx(iz,iy,ix,v)] 
                                     #if DIMENSIONS == 3
                                     * uh2lob[basisSize-1][idxConv(iz,z)]
                                     #endif
                                     * uh2lob[basisSize-1][idxConv(iy,y)] * uh2lob[basisSize-1][idxConv(ix,x)];
              }
            }
          }
        }
      }
    }
  }
  
  return lob;
}

/**
 * Auxilliary function to findMinMax
 * Project to GaussLobatto and modify the min/max if required
 */
void compareWithADERDGSolutionAtGaussLobattoNodes(const double* const luh, const int numberOfVariables, const int basisSize, double* const min, double* const max) {

#if DIMENSIONS == 3
  const int basisSize3D = basisSize;
#else
  const int basisSize3D = 1;  
#endif
const int order = basisSize-1;

  idx4 idx(basisSize3D, basisSize, basisSize, numberOfVariables);
  idx2 idxConv(basisSize, basisSize);
  
  double lobValue;
  int x,y,z,v,ix,iy,iz;
  
  for(z=0; z<basisSize3D; z++) {
    for(y=0; y<basisSize; y++) {
      for(x=0; x<basisSize; x++) {
        for(v=0; v<numberOfVariables; v++) {
          lobValue = 0;
          for(iz=0; iz<basisSize3D; iz++) {
            for(iy=0; iy<basisSize; iy++) {
              for(ix=0; ix<basisSize;ix++) {
                lobValue += luh[idx(iz,iy,ix,v)] 
                            #if DIMENSIONS == 3
                            * uh2lob[order][idxConv(iz,z)]
                            #endif
                            * uh2lob[order][idxConv(iy,y)] * uh2lob[order][idxConv(ix,x)];
              }
            }
          }
          min[v] = std::min( min[v], lobValue );
          max[v] = std::max( max[v], lobValue );
        }
      }
    }
  }

}

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel
