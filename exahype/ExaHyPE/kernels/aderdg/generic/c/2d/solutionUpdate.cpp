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

#include "../../../../GaussLegendreQuadrature.h"
#include "../../../../KernelUtils.h"

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

#if DIMENSIONS == 2

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int numberOfParameters,
                    const int basisSize) {
  const int order = basisSize - 1;

  idx3 idx(basisSize, basisSize, numberOfVariables);
  for (int j = 0; j < basisSize; j++) {
    for (int k = 0; k < basisSize; k++) {
      const double weight = kernels::gaussLegendreWeights[order][j] *
                            kernels::gaussLegendreWeights[order][k];
      const double updateSize = dt / weight;

      for (int l = 0; l < numberOfVariables - numberOfParameters; l++) {
        luh[idx(j, k, l)] += lduh[idx(j, k, l)] * updateSize;
      }
    }
  }
}

#endif  // DIMENSIONS == 2

}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
