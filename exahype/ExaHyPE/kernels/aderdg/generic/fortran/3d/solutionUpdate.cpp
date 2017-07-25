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

#include "string.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGMatrices.h"

using std::endl;
using std::cout;

extern "C" {
void elementupdate_(double *luh, double *lduh, double *dt);
}
void kernels::aderdg::generic::fortran::solutionUpdate(
    double *luh, const double *const lduh, const double dt,
    const int numberOfVariables, const int numberOfParameters, 
    const int basisSize) {

  double *lduhFortran =
      new double[numberOfVariables * basisSize * basisSize * basisSize];
  for (int i = 0; i < numberOfVariables * basisSize * basisSize * basisSize;
       i++) {
    lduhFortran[i] = lduh[i];
  }

  double *dtTemp = new double[1];
  dtTemp[0] = dt;

  elementupdate_(luh, lduhFortran, dtTemp);

  delete[] lduhFortran;
  delete dtTemp;
}
