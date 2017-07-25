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

#if DIMENSIONS == 3

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables,
                             const int numberOfParameters,
                             const int basisSize) {
  const int order = basisSize - 1;
  const int basisSize2 = basisSize * basisSize;
  const int basisSize3 = basisSize2 * basisSize;

  // Initialize the update DOF
  std::memset(lduh, 0, basisSize3 * numberOfVariables * sizeof(double));

  // x-direction
  {
    idx4 idx(basisSize, basisSize, basisSize, numberOfVariables);
    const int x_offset = 0 * basisSize3 * numberOfVariables;
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        const double weight = kernels::gaussLegendreWeights[order][i] *
                              kernels::gaussLegendreWeights[order][j];
        const double updateSize = weight / dx[0];

        // Fortran: lduh(l, k, j, i) += us * lFhi_x(l, m, j, i) * Kxi(m, k)
        // Matrix product: (l, m) * (m, k) = (l, k)
        for (int k = 0; k < basisSize; k++) {
          for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
            for (int m = 0; m < basisSize; m++) {
              lduh[idx(i, j, k, l)] += kernels::Kxi[order][k][m] *
                                       lFhi[x_offset + idx(i, j, m, l)] *
                                       updateSize;
            }
          }
        }
      }
    }
  }

  // y-direction
  {
    idx4 idx(basisSize, basisSize, basisSize, numberOfVariables);
    const int y_offset = 1 * basisSize3 * numberOfVariables;
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        const double weight = kernels::gaussLegendreWeights[order][i] *
                              kernels::gaussLegendreWeights[order][j];
        const double updateSize = weight / dx[1];

        // Fortran: lduh(l, j, k, i) += us * lFhi_y(l,m,j,i) * Kxi(m, k)
        // Matrix product: (l, m) * (m, k) = (l, k)
        for (int k = 0; k < basisSize; k++) {
          for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
            for (int m = 0; m < basisSize; m++) {
              lduh[idx(i, k, j, l)] += kernels::Kxi[order][k][m] *
                                       lFhi[y_offset + idx(i, j, m, l)] *
                                       updateSize;
            }
          }
        }
      }
    }
  }

  // z-direction
  {
    idx4 idx(basisSize, basisSize, basisSize, numberOfVariables);
    const int z_offset = 2 * basisSize3 * numberOfVariables;
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        const double weight = kernels::gaussLegendreWeights[order][i] *
                              kernels::gaussLegendreWeights[order][j];
        const double updateSize = weight / dx[2];

        // Fortran: lduh(l, j, i, k) += us * lFhi_z(l, m, j, i) * Kxi(m, k)
        // Matrix product (l, m) * (m, k) = (l, k)
        for (int k = 0; k < basisSize; k++) {
          for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
            for (int m = 0; m < basisSize; m++) {
              lduh[idx(k, i, j, l)] += kernels::Kxi[order][k][m] *
                                       lFhi[z_offset + idx(i, j, m, l)] *
                                       updateSize;
            }
          }
        }
      }
    }
  }

  // source
  {
    idx4 idx(basisSize, basisSize, basisSize, numberOfVariables);
    const int s_offset = 3 * basisSize3 * numberOfVariables;
    for (int i = 0; i < basisSize; i++) {
      for (int j = 0; j < basisSize; j++) {
        for (int k = 0; k < basisSize; k++) {
          const double weight = kernels::gaussLegendreWeights[order][i] *
                                kernels::gaussLegendreWeights[order][j] *
                                kernels::gaussLegendreWeights[order][k];

          // Fortran: lduh(:,k,j,i) += w * lShi(:,k,j,i)

          // TODO(guera): numberOfVariables - numberOfParameters
          for (int l = 0; l < numberOfVariables; l++) {
            lduh[idx(i, j, k, l)] += weight * lFhi[s_offset + idx(i, j, k, l)];
          }
        }
      }
    }
  }
}

#endif  // DIMENSIONS == 3

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
