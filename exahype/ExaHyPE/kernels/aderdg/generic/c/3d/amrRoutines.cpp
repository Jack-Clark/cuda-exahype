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
 
#include "kernels/aderdg/generic/Kernels.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

#include "string.h"

#if DIMENSIONS == 3
void singleLevelFaceUnknownsProlongation(
    double* lQhbndFine,
    const double* lQhbndCoarse,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m2 = 0; m2 < basisSize; ++m2) {
    for (int m1 = 0; m1 < basisSize; ++m1) {
      for (int ivar = 0; ivar < numberOfVariables; ++ivar) {
        const int mNodeIndex     = basisSize*m2 + m1;
        const int mDofStartIndex = mNodeIndex * numberOfVariables;

        for (int n2 = 0; n2 < basisSize; ++n2) {
          for (int n1 = 0; n1 < basisSize; ++n1) {
            const int nNodeIndex     = basisSize*n2 + n1;
            const int nDofStartIndex = nNodeIndex * numberOfVariables;

            lQhbndFine[mDofStartIndex+ivar] +=
                lQhbndCoarse[nDofStartIndex + ivar] *
                kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][n1][m1] *
                kernels::fineGridProjector1d[basisSize-1][subfaceIndex[1]][n2][m2];
          }
        }
      }
    }
  }
}

void kernels::aderdg::generic::c::faceUnknownsProlongation(double* lQhbndFine,
                                                           double* lFhbndFine,
                                                           const double* lQhbndCoarse,
                                                           const double* lFhbndCoarse,
                                                           const int coarseGridLevel,
                                                           const int fineGridLevel,
                                                           const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                                                           const int numberOfVariables, const int basisSize){
  const int levelDelta = fineGridLevel - coarseGridLevel;
  const int basisSize2 = basisSize*basisSize;

  double * lQhbndFineTemp = new double[basisSize2*numberOfVariables];
  double * lFhbndFineTemp = new double[basisSize2*numberOfVariables];

  double * pointerQhbnd1 = 0;
  double * pointerFhbnd1 = 0;

  double * pointerQhbnd2 = 0;
  double * pointerFhbnd2 = 0;

  // This ensures that the pointerQhbnd1 
  // of the last iteration points to lQhbndFine.
  // The same is done for pointerFhbnd1.
  if (levelDelta % 2 == 0) {
    pointerQhbnd1 = lQhbndFineTemp;
    pointerFhbnd1 = lFhbndFineTemp;
  } else {
    pointerQhbnd1 = lQhbndFine;
    pointerFhbnd1 = lFhbndFine;
  }

  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexPrevious (subfaceIndex);
  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexCurrent;

  tarch::la::Vector<DIMENSIONS-1,int> subintervalIndex;
  // This loop decodes the elements of subfaceIndex into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // 
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level prolongation.
  for (int l = 1; l < levelDelta+1; ++l) {
    const int significance = tarch::la::aPowI(levelDelta-l,3);
    subfaceIndexCurrent[0] = subfaceIndexPrevious[0] % significance;
    subfaceIndexCurrent[1] = subfaceIndexPrevious[1] % significance;
    subintervalIndex[0]    = (subfaceIndexPrevious[0] - subfaceIndexCurrent[0])/significance;
    subintervalIndex[1]    = (subfaceIndexPrevious[1] - subfaceIndexCurrent[1])/significance;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);

    // Zero the values of the first pointer.
    memset(pointerQhbnd1, 0, basisSize2*numberOfVariables*sizeof(double));
    memset(pointerFhbnd1, 0, basisSize2*numberOfVariables*sizeof(double));

    // Apply the single level prolongation operator.
    // Use the coarse level unknowns as input in the first iteration.
    if (l==1) {
      singleLevelFaceUnknownsProlongation(
          pointerQhbnd1,
          lQhbndCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsProlongation(
          pointerFhbnd1,
          lFhbndCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);

    } else {
      singleLevelFaceUnknownsProlongation(
          pointerQhbnd1,
          pointerQhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsProlongation(
          pointerFhbnd1,
          pointerFhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    subfaceIndexPrevious = subfaceIndexCurrent;

    pointerQhbnd2 = pointerQhbnd1;
    pointerFhbnd2 = pointerFhbnd1;
    
    // Toggle the addresses of the pointers.
    if (pointerQhbnd1 == lQhbndFineTemp) {
      pointerQhbnd1 = lQhbndFine;
      pointerFhbnd1 = lFhbndFine;
    } else {
      pointerQhbnd1 = lQhbndFineTemp;
      pointerFhbnd1 = lFhbndFineTemp;
    }
  }

  // Clean up.
  delete [] lQhbndFineTemp;
  delete [] lFhbndFineTemp;
}

void singleLevelFaceUnknownsRestriction(
    double* lQhbndCoarse,
    const double* lQhbndFine,
    const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m2 = 0; m2 < basisSize; ++m2) {
    for (int m1 = 0; m1 < basisSize; ++m1) {
      for (int ivar = 0; ivar < numberOfVariables; ++ivar) {
        const int mNodeIndex     = basisSize*m2 + m1;
        const int mDofStartIndex = mNodeIndex * numberOfVariables;

        for (int n2 = 0; n2 < basisSize; ++n2) {
          for (int n1 = 0; n1 < basisSize; ++n1) {
            const int nNodeIndex     = basisSize*n2 + n1;
            const int nDofStartIndex = nNodeIndex * numberOfVariables;
            
            lQhbndCoarse[mDofStartIndex+ivar] +=
                             kernels::gaussLegendreWeights[basisSize-1][n1] *
                             kernels::gaussLegendreWeights[basisSize-1][n2] *
                             kernels::fineGridProjector1d[basisSize-1][subfaceIndex[0]][m1][n1] *
                             kernels::fineGridProjector1d[basisSize-1][subfaceIndex[1]][m2][n2] *
                             lQhbndFine[nDofStartIndex + ivar] /
                             kernels::gaussLegendreWeights[basisSize-1][m1] /
                             kernels::gaussLegendreWeights[basisSize-1][m2] / 9.0;
          }
        }
      }
    }
  }
}

void accumulate(
    double * inOutArray,
    const double * inArray,
    const int length) {
  for (int i = 0; i < length; ++i) {
    inOutArray[i] += inArray[i];
  }
}

void kernels::aderdg::generic::c::faceUnknownsRestriction(double* lQhbndCoarse,
                                                          double* lFhbndCoarse,
                                                          const double* lQhbndFine,
                                                          const double* lFhbndFine,
                                                          const int coarseGridLevel,
                                                          const int fineGridLevel,
                                                          const tarch::la::Vector<DIMENSIONS-1, int>& subfaceIndex,
                                                          const int numberOfVariables, const int basisSize){
  const int levelDelta     = fineGridLevel - coarseGridLevel;
  const int basisSize2     = basisSize*basisSize;

  double * lQhbndCoarseTemp1 = new double[basisSize2*numberOfVariables];
  double * lFhbndCoarseTemp1 = new double[basisSize2*numberOfVariables];
  double * lQhbndCoarseTemp2 = new double[basisSize2*numberOfVariables];
  double * lFhbndCoarseTemp2 = new double[basisSize2*numberOfVariables];

  double * pointerQhbnd1 = 0;
  double * pointerQhbnd2 = 0;
  double * pointerFhbnd1 = 0;
  double * pointerFhbnd2 = 0;

  pointerQhbnd1 = lQhbndCoarseTemp1;
  pointerFhbnd1 = lFhbndCoarseTemp1;

  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndexCurrent(subfaceIndex);

  tarch::la::Vector<DIMENSIONS-1,int> subintervalIndex;
  
  // This loop decodes the indices of subfaceIndex into a tertiary basis
  // starting with the lowest significance 3^0 (in contrast to the prolongation loop).
  //
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level restriction.
  for (int l = 1; l < levelDelta+1; ++l) {
    subintervalIndex[0]    = subfaceIndexCurrent[0] % 3;
    subintervalIndex[1]    = subfaceIndexCurrent[1] % 3;
    subfaceIndexCurrent[0] = (subfaceIndexCurrent[0] - subintervalIndex[0])/3;
    subfaceIndexCurrent[1] = (subfaceIndexCurrent[1] - subintervalIndex[1])/3;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);

    // Zero the values of of the first pair of pointers.
    memset(pointerQhbnd1, 0, basisSize2*numberOfVariables*sizeof(double));
    memset(pointerFhbnd1, 0, basisSize2*numberOfVariables*sizeof(double));

    // Apply the single level restriction operator.
    // Use the fine level unknowns as input in the first iteration.
    if (l==1) {
      singleLevelFaceUnknownsRestriction(
          pointerQhbnd1,
          lQhbndFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsRestriction(
          pointerFhbnd1,
          lFhbndFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelFaceUnknownsRestriction(
          pointerQhbnd1,
          pointerQhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);

      singleLevelFaceUnknownsRestriction(
          pointerFhbnd1,
          pointerFhbnd2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    pointerQhbnd2 = pointerQhbnd1;
    pointerFhbnd2 = pointerFhbnd1;
    // Toggle pointer pairs.
    if (pointerQhbnd1 == lQhbndCoarseTemp1) {
      pointerQhbnd1 = lQhbndCoarseTemp2;
      pointerFhbnd1 = lFhbndCoarseTemp2;
    } else {
      pointerQhbnd1 = lQhbndCoarseTemp1;
      pointerFhbnd1 = lFhbndCoarseTemp1;
    }
  }

  // Add restricted fine level unknowns to coarse level unknowns.
  accumulate(lQhbndCoarse, pointerQhbnd2, basisSize2*numberOfVariables);
  accumulate(lFhbndCoarse, pointerFhbnd2, basisSize2*numberOfVariables);

  // Clean up.
  delete [] lQhbndCoarseTemp1;
  delete [] lFhbndCoarseTemp1;
  delete [] lQhbndCoarseTemp2;
  delete [] lFhbndCoarseTemp2;
}

void singleLevelVolumeUnknownsProlongation(
    double* luhFine,
    const double* luhCoarse,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m3 = 0; m3 < basisSize; ++m3) {
    for (int m2 = 0; m2 < basisSize; ++m2) {
      for (int m1 = 0; m1 < basisSize; ++m1) {
        for (int ivar = 0; ivar < numberOfVariables; ivar++) {
          const int mNodeIndex     = basisSize*basisSize*m3 + basisSize*m2 + m1;
          const int mDofStartIndex = mNodeIndex * numberOfVariables;

          for (int n3 = 0; n3 < basisSize; ++n3) {
            for (int n2 = 0; n2 < basisSize; ++n2) {
              for (int n1 = 0; n1 < basisSize; ++n1) {
                const int nNodeIndex     = basisSize*basisSize*n3 + basisSize*n2 + n1;
                const int nDofStartIndex = nNodeIndex * numberOfVariables;

                luhFine[mDofStartIndex+ivar]
                        += luhCoarse[nDofStartIndex + ivar] *
                        kernels::fineGridProjector1d[basisSize-1][subcellIndex[0]][n1][m1] *
                        kernels::fineGridProjector1d[basisSize-1][subcellIndex[1]][n2][m2] *
                        kernels::fineGridProjector1d[basisSize-1][subcellIndex[2]][n3][m3];
              }
            }
          }
        }
      }
    }
  }
}

void kernels::aderdg::generic::c::volumeUnknownsProlongation(double* luhFine,
                                                             const double* luhCoarse,
                                                             const int coarseGridLevel,
                                                             const int fineGridLevel,
                                                             const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
                                                             const int numberOfVariables, const int basisSize){
  const int levelDelta = fineGridLevel - coarseGridLevel;

  double * luhFineTemp = new double[basisSize*basisSize*basisSize*numberOfVariables];

  double * pointerUh1 = 0;
  double * pointerUh2 = 0;

  // This ensures that the first pointer 
  // points to luhFine in the last iteration
  // of the following loop.
  if (levelDelta % 2 == 0) {
    pointerUh1 = luhFineTemp;
  } else {
    pointerUh1 = luhFine;
  }

  tarch::la::Vector<DIMENSIONS,int> subcellIndexPrevious (subcellIndex);
  tarch::la::Vector<DIMENSIONS,int> subcellIndexCurrent;

  tarch::la::Vector<DIMENSIONS,int> subintervalIndex;
  // This loop step by step decodes the elements of subcellIndex into a tertiary basis
  // starting with the highest significance 3^(levelDelta-1).
  // 
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level prolongation.
  for (int l = 1; l < levelDelta+1; ++l) {
    const int significance = tarch::la::aPowI(levelDelta-l,3);
    subcellIndexCurrent[0] = subcellIndexPrevious[0] % significance;
    subcellIndexCurrent[1] = subcellIndexPrevious[1] % significance;
    subcellIndexCurrent[2] = subcellIndexPrevious[2] % significance;
    subintervalIndex[0]    = (subcellIndexPrevious[0] - subcellIndexCurrent[0])/significance;
    subintervalIndex[1]    = (subcellIndexPrevious[1] - subcellIndexCurrent[1])/significance;
    subintervalIndex[2]    = (subcellIndexPrevious[2] - subcellIndexCurrent[2])/significance;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);
    assertion(subintervalIndex[2] < 3);

    // Zero the values of the first pointer.
    memset(pointerUh1, 0, basisSize*basisSize*basisSize*numberOfVariables*sizeof(double));

    // Apply the single level prolongation operator.
    // Use the coarse level unknowns as input in the first iteration.
    if (l==1) {
      singleLevelVolumeUnknownsProlongation(
          pointerUh1,
          luhCoarse,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelVolumeUnknownsProlongation(
          pointerUh1,
          pointerUh2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    subcellIndexPrevious = subcellIndexCurrent;

    pointerUh2 = pointerUh1;

    // Toggle pointers.
    if (pointerUh1 == luhFineTemp) {
      pointerUh1 = luhFine;
    } else {
      pointerUh1 = luhFineTemp;
    }
  }

  // Clean up.
  delete [] luhFineTemp;
}

void singleLevelVolumeUnknownsRestriction(
    double* luhCoarse,
    const double* luhFine,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize) {
  for (int m3 = 0; m3 < basisSize; ++m3) {
    for (int m2 = 0; m2 < basisSize; ++m2) {
      for (int m1 = 0; m1 < basisSize; ++m1) {
        for (int ivar = 0; ivar < numberOfVariables; ivar++) {
          const int mNodeIndex     = basisSize*basisSize*m3 + basisSize*m2 + m1;
          const int mDofStartIndex = mNodeIndex * numberOfVariables;

          for (int n3 = 0; n3 < basisSize; ++n3) {
            for (int n2 = 0; n2 < basisSize; ++n2) {
              for (int n1 = 0; n1 < basisSize; ++n1) {
                const int nNodeIndex     = basisSize*basisSize*n3 + basisSize*n2 + n1;
                const int nDofStartIndex = nNodeIndex * numberOfVariables;

                luhCoarse[mDofStartIndex+ivar] +=
                    kernels::gaussLegendreWeights[basisSize-1][n1] *
                    kernels::gaussLegendreWeights[basisSize-1][n2] *
                    kernels::gaussLegendreWeights[basisSize-1][n3] *
                    kernels::fineGridProjector1d[basisSize-1][subcellIndex[0]][m1][n1] *
                    kernels::fineGridProjector1d[basisSize-1][subcellIndex[1]][m2][n2] *
                    kernels::fineGridProjector1d[basisSize-1][subcellIndex[2]][m3][n3] *
                    luhFine[nDofStartIndex + ivar] /
                    kernels::gaussLegendreWeights[basisSize-1][m1] /
                    kernels::gaussLegendreWeights[basisSize-1][m2] /
                    kernels::gaussLegendreWeights[basisSize-1][m3] / 27.0;
              }
            }
          }
        }
      }
    }
  }
}

void kernels::aderdg::generic::c::volumeUnknownsRestriction(
    double* luhCoarse,
    const double* luhFine,
    const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize){
  const int levelDelta     = fineGridLevel - coarseGridLevel;

  double * luhCoarseTemp1 = new double[basisSize*basisSize*basisSize*numberOfVariables];
  double * luhCoarseTemp2 = new double[basisSize*basisSize*basisSize*numberOfVariables];

  double * pointerUh1 = 0;
  double * pointerUh2 = 0;

  pointerUh1 = luhCoarseTemp1;

  tarch::la::Vector<DIMENSIONS,int> subcellIndexCurrent(subcellIndex);

  tarch::la::Vector<DIMENSIONS,int> subintervalIndex;
  // This loop step by step decodes the elements of subcellIndex into a tertiary basis
  // starting with the lowest significance 3^0 (in contrast to the prolongation decoding).
  //
  // Per iteration, the digits corresponding to the current significances then determine
  // the subintervals for the single level restriction.
  for (int l = 1; l < levelDelta+1; ++l) {
    subintervalIndex[0]    = subcellIndexCurrent[0] % 3;
    subintervalIndex[1]    = subcellIndexCurrent[1] % 3;
    subintervalIndex[2]    = subcellIndexCurrent[2] % 3;
    subcellIndexCurrent[0] = (subcellIndexCurrent[0] - subintervalIndex[0])/3;
    subcellIndexCurrent[1] = (subcellIndexCurrent[1] - subintervalIndex[1])/3;
    subcellIndexCurrent[2] = (subcellIndexCurrent[2] - subintervalIndex[2])/3;
    assertion(subintervalIndex[0] < 3);
    assertion(subintervalIndex[1] < 3);
    assertion(subintervalIndex[2] < 3);

    // Zero the values of the first pointer.
    memset(pointerUh1, 0, basisSize*basisSize*basisSize*numberOfVariables*sizeof(double));

    // Apply the single level restriction operator.
    // Use the fine level unknowns as input in the first iteration.
    if (l==1) {
      singleLevelVolumeUnknownsRestriction(
          pointerUh1,
          luhFine,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    } else {
      singleLevelVolumeUnknownsRestriction(
          pointerUh1,
          pointerUh2,
          subintervalIndex,
          numberOfVariables,
          basisSize);
    }

    // Prepare next iteration.
    pointerUh2 = pointerUh1;

    // Toggle the addresses of the pointers.
    if (pointerUh1 == luhCoarseTemp1) {
      pointerUh1 = luhCoarseTemp2;
    } else {
      pointerUh1 = luhCoarseTemp1;
    }
  }

  // Add restricted fine level unknowns to coarse level unknowns.
  accumulate(luhCoarse,pointerUh2,basisSize*basisSize*basisSize*numberOfVariables);

  // Clean up.
  delete [] luhCoarseTemp1;
  delete [] luhCoarseTemp2;
}
#endif
