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

#include <tarch/la/Vector.h>

#include "../../../../KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 2

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

void surfaceIntegralLinear(double *lduh, const double *const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double> &dx,
                           const int numberOfVariables, const int basisSize) {
  const int order = basisSize - 1;

  idx3 idx_lduh(basisSize, basisSize, numberOfVariables);
  idx3 idx_lFbnd(2 * DIMENSIONS, basisSize, numberOfVariables);

  // x faces
  for (int j = 0; j < basisSize; j++) {
    const double weight = kernels::gaussLegendreWeights[order][j];
    const double updateSize = weight / dx[0];

    for (int k = 0; k < basisSize; k++) {
      // left flux minus right flux
      for (int l = 0; l < numberOfVariables; l++) {
        lduh[idx_lduh(j, k, l)] -=
            (lFbnd[idx_lFbnd(1, j, l)] * kernels::FRCoeff[order][k] +
             lFbnd[idx_lFbnd(0, j, l)] * kernels::FLCoeff[order][k]) *
            updateSize;
      }
    }
  }

  // y faces
  for (int j = 0; j < basisSize; j++) {
    for (int k = 0; k < basisSize; k++) {
      const double weight = kernels::gaussLegendreWeights[order][k];
      const double updateSize = weight / dx[1];

      // back flux minus front flux
      for (int l = 0; l < numberOfVariables; l++) {
        lduh[idx_lduh(j, k, l)] -=
            (lFbnd[idx_lFbnd(3, k, l)] * kernels::FRCoeff[order][j] +
             lFbnd[idx_lFbnd(2, k, l)] * kernels::FLCoeff[order][j]) *
            updateSize;
      }
    }
  }
}

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#endif  // DIMENSIONS == 2
