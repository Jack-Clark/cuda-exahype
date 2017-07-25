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

#include "exahype/tests/kernels/c/ElasticityKernelTest.h"

#include <algorithm>
#include <cstring>
#include <limits>
#include <numeric>

#include "../testdata/elasticity_testdata.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 2

namespace exahype {
namespace tests {
namespace c {

void ElasticityKernelTest::flux(const double *Q, double **F) {
  for (int i = 0; i < kNumberOfVariables; i++) {
    assertion2(std::isfinite(Q[i]), Q[i], i);
  }

  double *f = F[0];
  double *g = F[1];
  // double* h = F[2];

  std::fill(f, f + kNumberOfVariables, 0.0);
  std::fill(g, g + kNumberOfVariables, 0.0);
  // std::fill(h, h + kNumberOfVariables, 0.0);

  double lam = Q[9];          // par(2)
  double mu = Q[10];          // par(2)
  double irho = 1.0 / Q[11];  // 1.0 / par(3)

  f[1 - 1] = -(lam + 2 * mu) * Q[7 - 1];
  f[2 - 1] = -lam * Q[7 - 1];
  f[3 - 1] = -lam * Q[7 - 1];
  f[4 - 1] = -mu * Q[8 - 1];
  f[5 - 1] = 0.0;
  f[6 - 1] = -mu * Q[9 - 1];
  f[7 - 1] = -irho * Q[1 - 1];
  f[8 - 1] = -irho * Q[4 - 1];
  f[9 - 1] = -irho * Q[6 - 1];

  g[1 - 1] = -lam * Q[8 - 1];
  g[2 - 1] = -(lam + 2 * mu) * Q[8 - 1];
  g[3 - 1] = -lam * Q[8 - 1];
  g[4 - 1] = -mu * Q[7 - 1];
  g[5 - 1] = -mu * Q[9 - 1];
  g[6 - 1] = 0.0;
  g[7 - 1] = -irho * Q[4 - 1];
  g[8 - 1] = -irho * Q[2 - 1];
  g[9 - 1] = -irho * Q[5 - 1];

  //  h[1 - 1] = - lam*Q[9 - 1];
  //  h[2 - 1] = - lam*Q[9 - 1];
  //  h[3 - 1] = - (lam+2*mu)*Q[9 - 1];
  //  h[4 - 1] = 0.0;
  //  h[5 - 1] = - mu *Q[8 - 1];
  //  h[6 - 1] = - mu *Q[7 - 1];
  //  h[7 - 1] = - irho *Q[4 - 1];
  //  h[8 - 1] = - irho *Q[2 - 1];
  //  h[9 - 1] = - irho *Q[5 - 1];

  for (int i = 0; i < DIMENSIONS; i++) {
    for (int j = 0; j < kNumberOfVariables; j++) {
      assertion3(std::isfinite(F[i][j]), F[i][j], i, j);
    }
  }
}

void ElasticityKernelTest::source(const double *Q, double *S) {
  std::fill(S, S + kNumberOfVariables, 0.0);
}

void ElasticityKernelTest::eigenvalues(const double *const Q,
                                           const int normalNonZeroIndex,
                                           double *lambda) {
  std::fill(lambda, lambda + kNumberOfParameters, 0.0);

  std::fill(lambda, lambda + kNumberOfVariables, 0.0);

  double lam = Q[9];    // par(1)
  double mu = Q[10];    // par(2)
  double rho0 = Q[11];  // par(3)
  double cp = std::sqrt((lam + 2 * mu) / rho0);
  double cs = std::sqrt(mu / rho0);

  assert(std::isfinite(cp));
  assert(std::isfinite(cs));

  lambda[1 - 1] = -cp;
  lambda[2 - 1] = -cs;
  lambda[3 - 1] = -cs;
  lambda[4 - 1] = 0.0;
  lambda[5 - 1] = 0.0;
  lambda[6 - 1] = 0.0;
  lambda[7 - 1] = +cs;
  lambda[8 - 1] = +cs;
  lambda[9 - 1] = +cp;
}

void ElasticityKernelTest::ncp(const double *const Q,
                                   const double *const gradQ, double *BgradQ) {
  std::fill(BgradQ, BgradQ + kNumberOfVariables * DIMENSIONS, 0.0);

  double lam = Q[kNumberOfVariables - kNumberOfParameters];     // par(1)
  double mu = Q[kNumberOfVariables - kNumberOfParameters + 1];  // par(2)
  double irho =
      1.0 / Q[kNumberOfVariables - kNumberOfParameters + 2];  // 1.0 / par(3)

  assert(std::isfinite(irho));

  const double *gradQx = gradQ + 0 * kNumberOfVariables;
  const double *gradQy = gradQ + 1 * kNumberOfVariables;
  //  const double *gradQz = gradQ + 2 * kNumberOfVariables;

  double *BgradQx = BgradQ + 0 * kNumberOfVariables;
  double *BgradQy = BgradQ + 1 * kNumberOfVariables;
  //  double *BgradQz = BgradQ + 2 * kNumberOfVariables;

  BgradQx[1 - 1] = -(lam + 2 * mu) * gradQx[7 - 1];
  BgradQx[2 - 1] = -lam * gradQx[7 - 1];
  BgradQx[3 - 1] = -lam * gradQx[7 - 1];
  BgradQx[4 - 1] = -mu * gradQx[8 - 1];
  BgradQx[5 - 1] = 0.0;
  BgradQx[6 - 1] = -mu * gradQx[9 - 1];
  BgradQx[7 - 1] = -irho * gradQx[1 - 1];
  BgradQx[8 - 1] = -irho * gradQx[4 - 1];
  BgradQx[9 - 1] = -irho * gradQx[6 - 1];

  BgradQy[1 - 1] = -lam * gradQy[8 - 1];
  BgradQy[2 - 1] = -(lam + 2 * mu) * gradQy[8 - 1];
  BgradQy[3 - 1] = -lam * gradQy[8 - 1];
  BgradQy[4 - 1] = -mu * gradQy[7 - 1];
  BgradQy[5 - 1] = -mu * gradQy[9 - 1];
  BgradQy[6 - 1] = 0.0;
  BgradQy[7 - 1] = -irho * gradQy[4 - 1];
  BgradQy[8 - 1] = -irho * gradQy[2 - 1];
  BgradQy[9 - 1] = -irho * gradQy[5 - 1];

  //  BgradQy[1 - 1] = -lam * gradQz[9 - 1];
  //  BgradQy[2 - 1] = -lam * gradQz[9 - 1];
  //  BgradQy[3 - 1] = -(lam + 2 * mu) * gradQz[9 - 1];
  //  BgradQy[4 - 1] = 0.0;
  //  BgradQy[5 - 1] = -mu * gradQz[8 - 1];
  //  BgradQy[6 - 1] = -mu * gradQz[7 - 1];
  //  BgradQy[7 - 1] = -irho * gradQz[6 - 1];
  //  BgradQy[8 - 1] = -irho * gradQz[5 - 1];
  //  BgradQy[9 - 1] = -irho * gradQz[3 - 1];

}  // ncp

void ElasticityKernelTest::matrixb(const double *const Q,
                                       const int normalNonZero, double *Bn) {
  std::fill(Bn, Bn + kNumberOfVariables * kNumberOfVariables, 0.0);

  kernels::idx2 idx_Bn(kNumberOfVariables, kNumberOfVariables);

  double lam = Q[9];          // par(1)
  double mu = Q[10];          // par(2)
  double irho = 1.0 / Q[11];  // 1./par(3)

  assert(std::isfinite(irho));

  switch (normalNonZero) {
    case 0:
      Bn[idx_Bn(7 - 1, 1 - 1)] = -(lam + 2 * mu);
      Bn[idx_Bn(7 - 1, 2 - 1)] = -lam;
      Bn[idx_Bn(7 - 1, 3 - 1)] = -lam;
      Bn[idx_Bn(8 - 1, 4 - 1)] = -mu;
      Bn[idx_Bn(9 - 1, 6 - 1)] = -mu;
      Bn[idx_Bn(1 - 1, 7 - 1)] = -irho;
      Bn[idx_Bn(4 - 1, 8 - 1)] = -irho;
      Bn[idx_Bn(6 - 1, 9 - 1)] = -irho;
      break;
    case 1:
      Bn[idx_Bn(8 - 1, 1 - 1)] = -lam;
      Bn[idx_Bn(8 - 1, 2 - 1)] = -(lam + 2 * mu);
      Bn[idx_Bn(8 - 1, 3 - 1)] = -lam;
      Bn[idx_Bn(7 - 1, 4 - 1)] = -mu;
      Bn[idx_Bn(9 - 1, 5 - 1)] = -mu;
      Bn[idx_Bn(4 - 1, 7 - 1)] = -irho;
      Bn[idx_Bn(2 - 1, 8 - 1)] = -irho;
      Bn[idx_Bn(5 - 1, 9 - 1)] = -irho;
      break;
    //    case 2:
    //      Bn[idx_Bn(9 - 1, 1 - 1)] = -lam;
    //      Bn[idx_Bn(9 - 1, 2 - 1)] = -lam;
    //      Bn[idx_Bn(9 - 1, 3 - 1)] = -(lam + 2 * mu);
    //      Bn[idx_Bn(8 - 1, 5 - 1)] = -mu;
    //      Bn[idx_Bn(7 - 1, 6 - 1)] = -mu;
    //      Bn[idx_Bn(6 - 1, 7 - 1)] = -irho;
    //      Bn[idx_Bn(5 - 1, 8 - 1)] = -irho;
    //      Bn[idx_Bn(3 - 1, 9 - 1)] = -irho;
    //      break;
    default:
      assert(false);
      break;
  }
}  // matrixb

void ElasticityKernelTest::testRiemannSolverLinear() {
  logInfo("ElasticityKernelTest::testRiemannSolverLinear()",
          "Test Riemann solver linear, ORDER=4, DIM=2");

  double *qL = new double[kBasisSize * kNumberOfVariables];
  double *qR = new double[kBasisSize * kNumberOfVariables];

  kernels::idx2 idx_q(kBasisSize, kNumberOfVariables);
  kernels::idx2 idx_q_in(kBasisSize, kNumberOfVariables - kNumberOfParameters);
  kernels::idx2 idx_param_in(kBasisSize, kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      qL[idx_q(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::qL_IN[idx_q_in(i, j)];
      qR[idx_q(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::qR_IN[idx_q_in(i, j)];
    }

    for (int j = 0; j < kNumberOfParameters; j++) {
      qL[idx_q(i, j + kNumberOfVariables - kNumberOfParameters)] =
          exahype::tests::testdata::elasticity::testRiemannSolverLinear::
              paramL_IN[idx_param_in(i, j)];
      qR[idx_q(i, j + kNumberOfVariables - kNumberOfParameters)] =
          exahype::tests::testdata::elasticity::testRiemannSolverLinear::
              paramR_IN[idx_param_in(i, j)];
    }
  }

  double *FL = new double[kBasisSize * kNumberOfVariables];
  double *FR = new double[kBasisSize * kNumberOfVariables];
  std::fill(FL, FL + kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());
  std::fill(FR, FR + kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());

  kernels::idx2 idx_F(kBasisSize, kNumberOfVariables);
  kernels::idx2 idx_F_in(kBasisSize, kNumberOfVariables - kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      FL[idx_F(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::FL_IN[idx_F_in(i, j)];
      FR[idx_F(i, j)] = exahype::tests::testdata::elasticity::
          testRiemannSolverLinear::FR_IN[idx_F_in(i, j)];
    }
  }

  const double dt = 1.916666666666667E-004;

  // temp variables:
  double  *tempFaceUnknowns      = new double[kBasisSize*kNumberOfVariables]; // nDOF(1) * nVar
  double **tempStateSizedVectors = new double*[5];
  tempStateSizedVectors[0] = new double[kNumberOfVariables*kNumberOfVariables]; // nVar
  tempStateSizedVectors[1] = tempStateSizedVectors[0]+kNumberOfVariables;
  tempStateSizedVectors[2] = tempStateSizedVectors[1]+kNumberOfVariables;
  tempStateSizedVectors[3] = tempStateSizedVectors[2]+kNumberOfVariables;
  tempStateSizedVectors[4] = tempStateSizedVectors[3]+kNumberOfVariables;
  double **tempStateSizedSquareMatrices = new double*[3];
  tempStateSizedSquareMatrices[0] = new double[3*kNumberOfVariables*kNumberOfVariables]; // nVar*nVar
  tempStateSizedSquareMatrices[1] = tempStateSizedSquareMatrices[0]+kNumberOfVariables*kNumberOfVariables;
  tempStateSizedSquareMatrices[2] = tempStateSizedSquareMatrices[1]+kNumberOfVariables*kNumberOfVariables;

  // TODO(Dominic): Fix test
  kernels::aderdg::generic::c::riemannSolverLinear<ElasticityKernelTest>(*this,
      FL, FR, qL, qR,
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,
      dt, 1 /* normalNonZero */);

  kernels::idx2 idx_F_out(kBasisSize, kNumberOfVariables - kNumberOfParameters);

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[idx_F(i, j)], exahype::tests::testdata::elasticity::
                               testRiemannSolverLinear::FL_OUT[idx_F_out(i, j)],
          eps, idx_F_out(i, j));
    }
  }

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kNumberOfVariables - kNumberOfParameters; j++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[idx_F(i, j)], exahype::tests::testdata::elasticity::
                               testRiemannSolverLinear::FR_OUT[idx_F_out(i, j)],
          eps, idx_F_out(i, j));
    }
  }

  delete[] qL;
  delete[] qR;
  delete[] FL;
  delete[] FR;
  delete[] tempStateSizedVectors[0];
  delete[] tempStateSizedVectors;
  delete[] tempStateSizedSquareMatrices[0];
  delete[] tempStateSizedSquareMatrices;
  delete[] tempFaceUnknowns;
}

void ElasticityKernelTest::testSpaceTimePredictorLinear() {
  logInfo("ElasticityKernelTest::testSpaceTimePredictorLinear()",
          "Test SpaceTimePredictor linear, ORDER=4, DIM=2");

  double *luh = new double[kNumberOfVariables * kBasisSize * kBasisSize];
  kernels::idx3 idx_luh(kBasisSize, kBasisSize, kNumberOfVariables);
  kernels::idx3 idx_luh_IN(kBasisSize, kBasisSize,
                           kNumberOfVariables - kNumberOfParameters);
  kernels::idx3 idx_param_IN(kBasisSize, kBasisSize, kNumberOfParameters);

  // Assemble luh = concatenate dofs and parameters
  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        luh[idx_luh(i, j, k)] = exahype::tests::testdata::elasticity::
            testSpaceTimePredictorLinear::luh_IN[idx_luh_IN(i, j, k)];
      }
      for (int k = 0; k < kNumberOfParameters; k++) {
        luh[idx_luh(i, j, k + kNumberOfVariables - kNumberOfParameters)] =
            exahype::tests::testdata::elasticity::testSpaceTimePredictorLinear::
                param_IN[idx_param_IN(i, j, k)];
      }
    }
  }

  // TODO: Unused
  // Inputs:
  double** tempSpaceTimeUnknowns = new double*[1];
  tempSpaceTimeUnknowns[0] = new double[kNumberOfVariables*kBasisSize*kBasisSize*(kBasisSize+1)];  // lQi; nVar * nDOFx * nDOFy * (nDOFt+1); nDOF+1 only here

  double** tempSpaceTimeFluxUnknowns = new double*[2];
  tempSpaceTimeFluxUnknowns[0] = new double[(DIMENSIONS+1)*kNumberOfVariables*kBasisSize*kBasisSize*kBasisSize];  // lFi+source; (dim+1) * nVar * nDOFx * nDOFy * nDOFt
  tempSpaceTimeFluxUnknowns[1] = new double[DIMENSIONS*kNumberOfVariables*kBasisSize*kBasisSize*kBasisSize];  // gradQ; dim * nVar * nDOFx * nDOFy * nDOFt

  double* tempStateSizedVector = nullptr;

  // Outputs:
  double *tempUnknowns     = new double[kNumberOfVariables*kBasisSize*kBasisSize*kBasisSize];     // lQhi; nVar * nDOFx * nDOFz
  double *tempFluxUnknowns = new double[(DIMENSIONS+1)*kNumberOfVariables*kBasisSize*kBasisSize]; // lFh+source,nVar * nDOFx * nDOFy * (dim+1)
  double *lQbnd = new double[kNumberOfVariables * kBasisSize * 2 * DIMENSIONS];
  double *lFbnd = new double[kNumberOfVariables * kBasisSize * 2 * DIMENSIONS];

  kernels::idx3 idx_lQhi(kBasisSize, kBasisSize, kNumberOfVariables);
  kernels::idx4 idx_lFhi(DIMENSIONS, kBasisSize, kBasisSize,
                         kNumberOfVariables);
  kernels::idx3 idx_lQbnd(2 * DIMENSIONS, kBasisSize, kNumberOfVariables);
  kernels::idx3 idx_lFbnd(2 * DIMENSIONS, kBasisSize, kNumberOfVariables);

  const tarch::la::Vector<DIMENSIONS, double> dx(38.4615384615385,
                                                 35.7142857142857);
  const double dt = 0.813172798364530;

  // Execute kernel
  // TODO(Dominic): Fix test
  kernels::aderdg::generic::c::spaceTimePredictorLinear<ElasticityKernelTest>(
      *this,
      lQbnd, lFbnd,
      tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
      tempUnknowns,tempFluxUnknowns,
      tempStateSizedVector,
      luh, dx, dt, nullptr);

  // Check result
  kernels::idx3 idx_lQhi_OUT(kBasisSize, kBasisSize,
                             kNumberOfVariables - kNumberOfParameters,
                             __LINE__);
  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        validateNumericalEqualsWithEpsWithParams1(
            tempUnknowns[idx_lQhi(i, j, k)],
            exahype::tests::testdata::elasticity::testSpaceTimePredictorLinear::
                lQhi_OUT[idx_lQhi_OUT(i, j, k)],
            eps, idx_lQhi_OUT(i, j, k));
      }
    }
  }

  kernels::idx4 idx_lFhi_OUT(DIMENSIONS, kBasisSize, kBasisSize,
                             kNumberOfVariables - kNumberOfParameters,
                             __LINE__);
  for (int i = 0; i < DIMENSIONS; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kBasisSize; k++) {
        for (int l = 0; l < kNumberOfVariables - kNumberOfParameters; l++) {
          validateNumericalEqualsWithEpsWithParams1(
              tempFluxUnknowns[idx_lFhi(i, j, k, l)],
              exahype::tests::testdata::elasticity::
                  testSpaceTimePredictorLinear::lFhi_OUT[idx_lFhi_OUT(i, j, k,
                                                                      l)],
              eps, idx_lFhi_OUT(i, j, k, l));
        }
      }
    }
  }

  kernels::idx3 idx_lQbnd_OUT(2 * DIMENSIONS, kBasisSize,
                              kNumberOfVariables - kNumberOfParameters,
                              __LINE__);
  for (int i = 0; i < 2 * DIMENSIONS; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        validateNumericalEqualsWithEpsWithParams1(
            lQbnd[idx_lQbnd(i, j, k)],
            exahype::tests::testdata::elasticity::testSpaceTimePredictorLinear::
                lQbnd_OUT[idx_lQbnd_OUT(i, j, k)],
            eps, idx_lQbnd_OUT(i, j, k));
      }
    }
  }

  kernels::idx3 idx_lFbnd_OUT(2 * DIMENSIONS, kBasisSize,
                              kNumberOfVariables - kNumberOfParameters,
                              __LINE__);
  for (int i = 0; i < 2 * DIMENSIONS; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        validateNumericalEqualsWithEpsWithParams1(
            lFbnd[idx_lFbnd(i, j, k)],
            exahype::tests::testdata::elasticity::testSpaceTimePredictorLinear::
                lFbnd_OUT[idx_lFbnd_OUT(i, j, k)],
            eps, idx_lFbnd_OUT(i, j, k));
      }
    }
  }

  delete[] luh;
  delete[] lQbnd;
  delete[] lFbnd;

  delete[] tempSpaceTimeUnknowns[0];
  delete[] tempSpaceTimeUnknowns;

  delete[] tempSpaceTimeFluxUnknowns[0];
  delete[] tempSpaceTimeFluxUnknowns[1];
  delete[] tempSpaceTimeFluxUnknowns;

  delete[] tempUnknowns;
  delete[] tempFluxUnknowns;
}

void ElasticityKernelTest::testVolumeIntegralLinear() {
  logInfo("ElasticityKernelTest::testVolumeIntegralLinear()",
          "Test VolumeIntegral linear, ORDER=4, DIM=2");

  double *lFhi = new double[(DIMENSIONS + 1) * kBasisSize * kBasisSize *
                            kNumberOfVariables];
  std::fill(
      lFhi,
      lFhi + (DIMENSIONS + 1) * kBasisSize * kBasisSize * kNumberOfVariables,
      std::numeric_limits<double>::quiet_NaN());
  kernels::idx4 idx_lFhi(DIMENSIONS + 1, kBasisSize, kBasisSize,
                         kNumberOfVariables);
  kernels::idx4 idx_lFhi_IN(DIMENSIONS, kBasisSize, kBasisSize,
                            kNumberOfVariables - kNumberOfParameters);

  for (int i = 0; i < DIMENSIONS; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kBasisSize; k++) {
        for (int l = 0; l < kNumberOfVariables - kNumberOfParameters; l++) {
          lFhi[idx_lFhi(i, j, k, l)] = exahype::tests::testdata::elasticity::
              testVolumeIntegralLinear::lFhi_IN[idx_lFhi_IN(i, j, k, l)];
        }
      }
    }
  }

  double *lduh = new double[kBasisSize * kBasisSize * kNumberOfVariables];
  kernels::idx3 idx_lduh(kBasisSize, kBasisSize, kNumberOfVariables);
  kernels::idx3 idx_lduh_OUT(kBasisSize, kBasisSize,
                             kNumberOfVariables - kNumberOfParameters);

  const tarch::la::Vector<DIMENSIONS, double> dx(38.4615384615385,
                                                 35.7142857142857);

  // Execute kernel
  kernels::aderdg::generic::c::volumeIntegralLinear(
      lduh, lFhi, dx, kNumberOfVariables, kNumberOfParameters, kBasisSize);

  // Compare
  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        validateNumericalEqualsWithEpsWithParams1(
            lduh[idx_lduh(i, j, k)],
            exahype::tests::testdata::elasticity::testVolumeIntegralLinear::
                lduh_OUT[idx_lduh_OUT(i, j, k)],
            eps, idx_lduh_OUT(i, j, k));
      }
    }
  }

  delete[] lFhi;
  delete[] lduh;
}

void ElasticityKernelTest::testSurfaceIntegralLinear() {
  logInfo("ElasticityKernelTest::testSurfaceIntegralLinear()",
          "Test SurfaceIntegral linear, ORDER=4, DIM=2");

  double *lFbnd = new double[2 * DIMENSIONS * kBasisSize * kNumberOfVariables];
  kernels::idx3 idx_lFbnd(2 * DIMENSIONS, kBasisSize, kNumberOfVariables);
  kernels::idx3 idx_lFbnd_IN(2 * DIMENSIONS, kBasisSize,
                             kNumberOfVariables - kNumberOfParameters);

  std::fill(lFbnd, lFbnd + 2 * DIMENSIONS * kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());

  for (int i = 0; i < 2 * DIMENSIONS; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        lFbnd[idx_lFbnd(i, j, k)] = exahype::tests::testdata::elasticity::
            testSurfaceIntegralLinear::lFbnd_IN[idx_lFbnd_IN(i, j, k)];
      }
    }
  }

  double *lduh = new double[kBasisSize * kBasisSize * kNumberOfVariables];
  kernels::idx3 idx_lduh(kBasisSize, kBasisSize, kNumberOfVariables);
  kernels::idx3 idx_lduh_INOUT(kBasisSize, kBasisSize,
                               kNumberOfVariables - kNumberOfParameters);

  std::fill(lduh, lduh + kBasisSize * kBasisSize * kNumberOfVariables,
            std::numeric_limits<double>::quiet_NaN());

  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        lduh[idx_lduh(i, j, k)] = exahype::tests::testdata::elasticity::
            testSurfaceIntegralLinear::lduh_IN[idx_lduh_INOUT(i, j, k)];
      }
    }
  }

  const tarch::la::Vector<DIMENSIONS, double> dx(38.4615384615385,
                                                 35.7142857142857);

  // Execute kernel
  kernels::aderdg::generic::c::surfaceIntegralLinear(
      lduh, lFbnd, dx, kNumberOfVariables, kBasisSize);

  // Compare
  for (int i = 0; i < kBasisSize; i++) {
    for (int j = 0; j < kBasisSize; j++) {
      for (int k = 0; k < kNumberOfVariables - kNumberOfParameters; k++) {
        validateNumericalEqualsWithEpsWithParams1(
            lduh[idx_lduh(i, j, k)],
            exahype::tests::testdata::elasticity::testSurfaceIntegralLinear::
                lduh_OUT[idx_lduh_INOUT(i, j, k)],
            eps, idx_lduh_INOUT(i, j, k));
      }
    }
  }

  delete[] lFbnd;
  delete[] lduh;
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==2
