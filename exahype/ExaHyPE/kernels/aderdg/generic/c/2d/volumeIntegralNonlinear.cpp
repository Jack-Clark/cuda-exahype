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

#include "../../Kernels.h"

#include <cstring>

#include <tarch/la/Vector.h>

#include "../../../../DGMatrices.h"
#include "../../../../GaussLegendreQuadrature.h"
#include "../../../../KernelUtils.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 2

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables,
                             const int numberOfParameters,
                             const int basisSize) {
  const int order = basisSize - 1;
  const int basisSize2 = basisSize * basisSize;

  // Initialize the update DOF
  std::memset(lduh, 0, basisSize2 * numberOfVariables * sizeof(double));

  // x-direction
  {
    idx3 idx(basisSize, basisSize, numberOfVariables);
    const int x_offset = 0 * basisSize2 * numberOfVariables;
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[0];

      // Fortran: lduh(l, k, j) += lFhi_x(l, m, j) * Kxi(m, k)
      // Matrix product: (l, m) * (m, k) = (l, k)
      for (int k = 0; k < basisSize; k++) {
        for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
          for (int m = 0; m < basisSize; m++) {
            lduh[idx(j, k, l)] += kernels::Kxi[order][k][m] *
                                  lFhi[x_offset + idx(j, m, l)] * updateSize;
          }
        }
      }
    }
  }

  // y-direction
  {
    idx3 idx(basisSize, basisSize, numberOfVariables);
    const int y_offset = 1 * basisSize2 * numberOfVariables;
    for (int j = 0; j < basisSize; j++) {
      const double weight = kernels::gaussLegendreWeights[order][j];
      const double updateSize = weight / dx[1];

      // Fortran: lduh(l, j, k) += lFhi_y(l, m, j) * Kxi(m, k)
      // Matrix product: (l, m) * (m, k) = (l, k)
      for (int k = 0; k < basisSize; k++) {
        for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
          for (int m = 0; m < basisSize; m++) {
            lduh[idx(k, j, l)] += kernels::Kxi[order][k][m] *
                                  lFhi[y_offset + idx(j, m, l)] * updateSize;
          }
        }
      }
    }
  }

  // source
  {
    idx3 idx(basisSize, basisSize, numberOfVariables);
    const int s_offset = 2 * basisSize2 * numberOfVariables;
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
        const double weight = kernels::gaussLegendreWeights[order][j] *
                              kernels::gaussLegendreWeights[order][k];

        // Fortran: lduh(:,k,j) += w * lShi(:,k,j)

        // TODO(guera): numberOfVariables - numberOfParameters
        for (int l = 0; l < numberOfVariables; l++) {
          lduh[idx(j, k, l)] += weight * lFhi[s_offset + idx(j, k, l)];
        }
      }
    }
  }
}

// TODO(Dominic): Commented out for now.
//void spaceTimeVolumeIntegralNonlinear(double* lduh, const double* const lFi,
//                             const tarch::la::Vector<DIMENSIONS, double>& dx,
//                             const int numberOfVariables,
//                             const int numberOfParameters,
//                             const int basisSize) {
//  const int order = basisSize - 1;
//  const int basisSize2 = basisSize * basisSize * basisSize;
//  const int basisSize3 = basisSize2 * basisSize;
//
//  // Initialize the update DOF
//  std::memset(lduh, 0, basisSize2 * numberOfVariables * sizeof(double));
//
//  for (int i = 0; i < basisSize; i++) { // loop over time unknowns
//    // x-direction
//    {
//      idx3 idx(basisSize, basisSize, numberOfVariables);
//      idx4 idx_lFix(basisSize,basisSize,basisSize,numberOfVariables);
//      const int x_offset = 0 * basisSize3 * numberOfVariables;
//      for (int j = 0; j < basisSize; j++) {
//        const double weight =
//            kernels::gaussLegendreWeights[order][i] *
//            kernels::gaussLegendreWeights[order][j];
//        const double updateSize = weight / dx[0];
//
//        // Fortran: lduh(l, k, j) += lFhi_x(l, m, j) * Kxi(m, k)
//        // Matrix product: (l, m) * (m, k) = (l, k)
//        for (int k = 0; k < basisSize; k++) {
//          for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
//            for (int m = 0; m < basisSize; m++) {
//              lduh[idx(j, k, l)] += kernels::Kxi[order][k][m] *
//                  lFi[x_offset + idx_lFix(i, j, m, l)] * updateSize;
//            }
//          }
//        }
//      }
//    }
//
//    // y-direction
//    {
//      idx3 idx(basisSize, basisSize, numberOfVariables);
//      idx4 idx_lFiy(basisSize,basisSize,basisSize,numberOfVariables);
//      const int y_offset = 1 * basisSize3 * numberOfVariables;
//      for (int j = 0; j < basisSize; j++) {
//        const double weight =
//            kernels::gaussLegendreWeights[order][i] *
//            kernels::gaussLegendreWeights[order][j];
//        const double updateSize = weight / dx[1];
//
//        // Fortran: lduh(l, j, k) += lFhi_y(l, m, j) * Kxi(m, k)
//        // Matrix product: (l, m) * (m, k) = (l, k)
//        for (int k = 0; k < basisSize; k++) {
//          for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
//            for (int m = 0; m < basisSize; m++) {
//              lduh[idx(k, j, l)] += kernels::Kxi[order][k][m] *
//                  lFi[y_offset + idx_lFiy(i, j, m, l)] * updateSize;
//            }
//          }
//        }
//      }
//    }
//
//    // source
//    {
//      idx3 idx(basisSize, basisSize, numberOfVariables);
//      idx4 idx_lS(basisSize,basisSize,basisSize,numberOfVariables);
//      const int s_offset = 2 * basisSize3 * numberOfVariables;
//      for (int j = 0; j < basisSize; j++) {
//        for (int k = 0; k < basisSize; k++) {
//          const double weight =
//              kernels::gaussLegendreWeights[order][i] *
//              kernels::gaussLegendreWeights[order][j] *
//              kernels::gaussLegendreWeights[order][k];
//
//          // Fortran: lduh(:,k,j) += w * lShi(:,k,j)
//
//          // TODO(guera): numberOfVariables - numberOfParameters
//          for (int l = 0; l < numberOfVariables; l++) {
//            lduh[idx(j, k, l)] += weight * lFi[s_offset + idx_lS(i, j, k, l)];
//          }
//        }
//      }
//    }
//  }
//}

#endif  // DIMENSIONS == 2

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
