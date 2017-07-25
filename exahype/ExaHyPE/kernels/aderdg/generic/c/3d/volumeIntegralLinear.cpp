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

#include "../../../../KernelUtils.h"
#include <cstring>

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 3

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables,
                          const int numberOfParameters, const int basisSize) {
  // for linear non-conservative PDE, the volume integral is trivial, since it
  // only involves the element mass matrix, which later will cancel

  idx5 idx_lFhi(DIMENSIONS, basisSize, basisSize, basisSize, numberOfVariables);
  idx4 idx_lduh(basisSize, basisSize, basisSize, numberOfVariables);

  const int basisSize2 = basisSize * basisSize;
  const int basisSize3 = basisSize2 * basisSize;
  const int order = basisSize - 1;

  std::memset(lduh, 0, basisSize3 * numberOfVariables * sizeof(double));

  for (int i = 0; i < basisSize; i++) {
    for (int j = 0; j < basisSize; j++) {
      for (int k = 0; k < basisSize; k++) {
        double weight = kernels::gaussLegendreWeights[order][i] *
                        kernels::gaussLegendreWeights[order][j] *
                        kernels::gaussLegendreWeights[order][k];

        // Fortran: lduh(:,k,j,i) = -SUM(lFhi(:,k,j,i,1:nDim), dim = 5) * weight
        // Avoid diffusion of parameters
        for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
	  double sum = 0.0;
          for (int m = 0; m < DIMENSIONS; m++) {
	    sum = sum  +   lFhi[idx_lFhi(m, i, j, k, l)];
            //lduh[idx_lduh(i, j, k, l)] -=
            //    weight * lFhi[idx_lFhi(m, i, j, k, l)];
          }
	  lduh[idx_lduh(i, j, k, l)] = -weight *sum;
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
