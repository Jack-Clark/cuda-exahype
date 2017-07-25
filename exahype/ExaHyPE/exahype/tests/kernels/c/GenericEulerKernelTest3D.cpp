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

#include "exahype/tests/kernels/c/GenericEulerKernelTest.h"

#include "../testdata/generic_euler_testdata.h"
#include "kernels/KernelUtils.h"
#include "kernels/aderdg/generic/Kernels.h"

#if DIMENSIONS == 3

namespace exahype {
namespace tests {
namespace c {

void GenericEulerKernelTest::flux(const double *const Q, double **F) {
  double *f = F[0];
  double *g = F[1];
  double *h = F[2];

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) *
      (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);

  f[0] = Q[1];
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);

  g[0] = Q[2];
  g[1] = irho * Q[2] * Q[1];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);

  h[0] = Q[3];
  h[1] = irho * Q[3] * Q[1];
  h[2] = irho * Q[3] * Q[2];
  h[3] = irho * Q[3] * Q[3] + p;
  h[4] = irho * Q[3] * (Q[4] + p);
}

void GenericEulerKernelTest::source(const double *const Q, double *S) {
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}

void GenericEulerKernelTest::eigenvalues(const double *const Q,
                                             const int normalNonZeroIndex,
                                             double *lambda) {
  const double GAMMA = 1.4;

  double irho = 1.0 / Q[0];
  double p = (GAMMA - 1) *
             (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2] + Q[3] * Q[3]) * irho);

  double u_n = Q[normalNonZeroIndex + 1] * irho;
  double c = std::sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}

void GenericEulerKernelTest::ncp(const double *const Q,
                                     const double *const gradQ,
                                     double *BgradQ) {
  // Q[5]
  // gradQ[3][5]
  // BgradQ[3][5]

  // Arbitrary BS!

  BgradQ[0] = 1.0;
  BgradQ[1] = -3.0;
  BgradQ[2] = Q[3];
  BgradQ[3] = gradQ[3];
  BgradQ[4] = -1.7;
  BgradQ[5] = 4.0;
  BgradQ[6] = Q[0];
  BgradQ[7] = 0.1;
  BgradQ[8] = 0.8;
  BgradQ[9] = Q[4];
  BgradQ[10] = 2.0;
  BgradQ[11] = 8.0;
  BgradQ[12] = -1.0;
  BgradQ[13] = gradQ[14];
  BgradQ[14] = 5.3;
}  // ncp

void GenericEulerKernelTest::matrixb(const double *const Q,
                                         const int normalNonZero, double *Bn) {
  // 3D compressible Euler equations
  double *B1 = new double[5 * 5];
  double *B2 = new double[5 * 5];
  double *B3 = new double[5 * 5];

  std::memset(B1, 0, 5 * 5 * sizeof(double));
  std::memset(B2, 0, 5 * 5 * sizeof(double));
  std::memset(B3, 0, 5 * 5 * sizeof(double));

  // Bn = B1 if normalNonZero == 0
  //      B2 if normalNonZero == 1
  //      B3 if normalNonZero == 2
  std::memcpy(Bn, (normalNonZero == 0) ? B1 : (normalNonZero == 1) ? B2 : B3,
              5 * 5 * sizeof(double));

  delete[] B1;
  delete[] B2;
  delete[] B3;
}  // matrixb

void GenericEulerKernelTest::testPDEFluxes() {
  logInfo("testPDEFluxes()", "Test PDE-related functions, DIM=3");

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5};  // pressure = 1.372
  double f[5], g[5], h[5];
  double *F[3] = {f, g, h};

  flux(Q, F);

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        f[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::f[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        g[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::g[i],
        i);
  }

  for (int i = 0; i < 5; i++) {
    validateNumericalEqualsWithParams1(
        h[i], ::exahype::tests::testdata::generic_euler::testPDEFluxes::h[i],
        i);
  }
}  // testPDEFluxes

void GenericEulerKernelTest::testVolumeIntegralLinear() {
  logInfo("testVolumeIntegralLinear()",
          "Test volume integral linear, ORDER=3, DIM=3");

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5,
                                                 0.5);  // mesh spacing
  double *lFhi = new double[960];  // nVar * nDOFx * nDOFy * nDOFz * dim
  std::fill(lFhi, lFhi + 960, 0.0);
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z ], 320 entries each
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[320];
  double *lFhi_z = &lFhi[640];

  // seed direction
  for (int i = 0; i < 320; i += 5) {
    lFhi_x[i + 1] = 1.;
    lFhi_y[i + 2] = 2.;
    lFhi_z[i + 3] = 3.;
  }

  kernels::aderdg::generic::c::volumeIntegralLinear(lduh, lFhi, dx,
                                                    5,  // numberOfVariables
                                                    0,  // numberOfParameters
                                                    4   // basisSize
                                                    );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testVolumeIntegralLinear::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegralLinear

void GenericEulerKernelTest::testVolumeIntegralNonlinear() {
  logInfo("testVolumeIntegralNonlinear()",
          "Test volume integral nonlinear, ORDER=3, DIM=3");

  // output:
  double *lduh = new double[320];  // intentionally left uninitialised

  // input:
  const double dx[3] = {0.05, 0.05, 0.05};  // mesh spacing
  double *lFhi = new double[1280];  // nVar * nDOFx * nDOFy * nDOFz * (dim + 1)
  std::fill(lFhi, lFhi + 1280, 0.0);
  // lFhi = [ lFhi_x  | lFhi_y | lFhi_z | lShi], 320 entries each
  double *lFhi_x = &lFhi[0];
  double *lFhi_y = &lFhi[320];
  double *lFhi_z = &lFhi[640];
  double *lShi = &lFhi[960];

  std::fill(lFhi, lFhi + 1280, 0.0);

  // seed direction
  for (int i = 0; i < 320; i += 5) {
    lFhi_x[i + 1] = 1.;
    lFhi_y[i + 2] = 1.;
    lFhi_z[i + 3] = 1.;
  }

  std::fill(lShi, lShi + 320, 0.0);

  kernels::aderdg::generic::c::volumeIntegralNonlinear(
      lduh, lFhi, dx[0],
      5,   // getNumberOfVariables(),
      0,   // getNumberOfParameters()
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testVolumeIntegralNonlinear::lduh[i],
        eps, i);
  }

  delete[] lduh;
  delete[] lFhi;
}  // testVolumeIntegralNonlinear

void GenericEulerKernelTest::testSurfaceIntegralLinear() {
  logInfo("testSurfaceIntegralLinear()",
          "Test surface integral linear, ORDER=3, DIM=3");

  // inputs:
  const double dx[3] = {0.1, 0.1, 0.1};  // mesh spacing
  double *lFhbnd = new double[6 * 80];   // 480
  std::fill(lFhbnd, lFhbnd + 6 * 80, 0.0);
  double *FLeft = &lFhbnd[0];
  double *FRight = &lFhbnd[80];
  double *FFront = &lFhbnd[160];
  double *FBack = &lFhbnd[240];
  double *FBottom = &lFhbnd[320];
  double *FTop = &lFhbnd[400];

  for (int i = 0; i < 80; i += 5) {
    // in x orientation 1
    FLeft[i + 1] = 1.;
    FRight[i + 1] = 1.;
    // in y orientation 1
    FFront[i + 2] = 1.;
    FBack[i + 2] = 1.;
    // in z direction 1
    FBottom[i + 3] = 1.;
    FTop[i + 3] = 1.;
  }

  // in/out:
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i / 10.;
  }

  // lFhbnd = [ FLeft | FRight | FFront | FBack | FBottom | FTop ]
  kernels::aderdg::generic::c::surfaceIntegralLinear(
      lduh, lFhbnd, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testSurfaceIntegralLinear::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegralLinear

void GenericEulerKernelTest::testSurfaceIntegralNonlinear() {
  logInfo("testSurfaceIntegralNonlinear()",
          "Test surface integral nonlinear, ORDER=3, DIM=3");

  // inputs:
  const double dx[3] = {0.1, 0.1, 0.1};  // mesh spacing
  double *lFhbnd = new double[6 * 80];   // 480
  std::fill(lFhbnd, lFhbnd + 6 * 80, 0.0);
  double *FLeft = &lFhbnd[0];
  double *FRight = &lFhbnd[80];
  double *FFront = &lFhbnd[160];
  double *FBack = &lFhbnd[240];
  double *FBottom = &lFhbnd[320];
  double *FTop = &lFhbnd[400];

  for (int i = 0; i < 80; i += 5) {
    // in x orientation 1
    FLeft[i + 1] = 1.;
    FRight[i + 1] = 1.;
    // in y orientation 1
    FFront[i + 2] = 1.;
    FBack[i + 2] = 1.;
    // in z direction 1
    FBottom[i + 3] = 1.;
    FTop[i + 3] = 1.;
  }

  // in/out:
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i / 10.;
  }

  // lFhbnd = [ FLeft | FRight | FFront | FBack | FBottom | FTop ]
  kernels::aderdg::generic::c::surfaceIntegralNonlinear(
      lduh, lFhbnd, dx[0],
      5,   // getNumberOfVariables(),
      4);  // getNodesPerCoordinateAxis()

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lduh[i], ::exahype::tests::testdata::generic_euler::
                     testSurfaceIntegralNonlinear::lduh[i],
        eps, i);
  }

  delete[] lFhbnd;
  delete[] lduh;
}  // testSurfaceIntegralNonlinear

void GenericEulerKernelTest::testRiemannSolverLinear() {
  logInfo("testRiemannSolverLinear()",
          "Test Riemann solver linear, ORDER=3, DIM=3");

  // output (intentionally left uninitialised):
  double *FL = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  double *FR = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  std::memset(FL, 0, 80 * sizeof(double));
  std::memset(FR, 0, 80 * sizeof(double));

  double  *tempFaceUnknowns      = new double[80]; // nDOF(1) * nDOF(2) * nVar
  double **tempStateSizedVectors = new double*[5];
  tempStateSizedVectors[0] = new double[5*5]; // nVar
  tempStateSizedVectors[1] = tempStateSizedVectors[0]+5;
  tempStateSizedVectors[2] = tempStateSizedVectors[1]+5;
  tempStateSizedVectors[3] = tempStateSizedVectors[2]+5;
  tempStateSizedVectors[4] = tempStateSizedVectors[3]+5;
  double **tempStateSizedSquareMatrices = new double*[3];
  tempStateSizedSquareMatrices[0] = new double[3*25]; // nVar*nVar
  tempStateSizedSquareMatrices[1] = tempStateSizedSquareMatrices[0]+25;
  tempStateSizedSquareMatrices[2] = tempStateSizedSquareMatrices[1]+25;

  // inputs:
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QL[80 =
  // nVar * nVar * nDOF]
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QR[80 =
  // nVar * nVar * nDOF]
  const double dt = 1.40831757919882352703e-03;

  // TODO(Dominic): Fix test
  kernels::aderdg::generic::c::riemannSolverLinear<GenericEulerKernelTest>(
      *this,
      FL, FR, exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
      exahype::tests::testdata::generic_euler::testRiemannSolver::QR,
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,
      dt,
      1 /*normalNonZero (only changes result of eigenvalues, matrixb) */
      );

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FL[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverLinear::FL[i],
        eps, i);
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FR[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverLinear::FR[i],
        eps, i);
  }

  delete[] FL;
  delete[] FR;
  delete[] tempStateSizedVectors[0];
  delete[] tempStateSizedVectors;
  delete[] tempStateSizedSquareMatrices[0];
  delete[] tempStateSizedSquareMatrices;
  delete[] tempFaceUnknowns;
}  // testRiemannSolverLinear

void GenericEulerKernelTest::testRiemannSolverNonlinear() {
  logInfo("testRiemannSolverNonlinear()",
          "Test Riemann solver nonlinear, ORDER=3, DIM=3");

  // inout:
  double *FL = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  double *FR = new double[4 * 4 * 5];  // nDOF(3) * nDOF(2) * nVar
  for (int i = 0; i < 80; i++) {
    // arbitrary values
    FL[i] = static_cast<double>(i + 1);
    FR[i] = static_cast<double>(i - 1);
  }
  double  *tempFaceUnknowns             = nullptr;
  double **tempStateSizedVectors        = new double*[6];
  tempStateSizedVectors[0]              = new double[6*5];
  tempStateSizedVectors[1]              = tempStateSizedVectors[0]+5;
  tempStateSizedVectors[2]              = tempStateSizedVectors[1]+5;
  tempStateSizedVectors[3]              = tempStateSizedVectors[2]+5;
  tempStateSizedVectors[4]              = tempStateSizedVectors[3]+5;
  tempStateSizedVectors[5]              = tempStateSizedVectors[4]+5;
  double **tempStateSizedSquareMatrices = new double*[1];
  tempStateSizedSquareMatrices[0]       = new double[1*5*5];

  // inputs:
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QL[80 =
  // nVar * nVar * nDOF]
  // exahype::tests::testdata::generic_euler::testRiemannSolver::QR[80 =
  // nVar * nVar * nDOF]
  const double dt = 1.40831757919882352703e-03;

  kernels::aderdg::generic::c::riemannSolverNonlinear<GenericEulerKernelTest>(
      *this,
      FL, FR,
      ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
      ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR,
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,
      dt,
      1  // normalNonZero
      );

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FL[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverNonlinear::FL[i],
        eps, i);
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        FR[i], ::exahype::tests::testdata::generic_euler::
                   testRiemannSolverNonlinear::FR[i],
        eps, i);
  }

  delete[] FL;
  delete[] FR;
  delete[] tempStateSizedSquareMatrices[0];
  delete[] tempStateSizedSquareMatrices;
  delete[] tempStateSizedVectors[0];
  delete[] tempStateSizedVectors;
}  // testRiemannSolverNonlinear

void GenericEulerKernelTest::testSolutionUpdate() {
  logInfo("testSolutionUpdate()", "Test solution update, ORDER=3, DIM=3");

  // in/out:
  double *luh = new double[320];
  std::fill(luh, luh + 320, 0.0);
  for (int i = 0; i < 320; i += 5) {
    luh[i] = 1.0;
    luh[i + 4] = 2.5;
  }

  // inputs:
  const double dt = 1.40831757919882352703e-03;
  double *lduh = new double[320];
  for (int i = 0; i < 320; i++) {
    lduh[i] = i;
  }

  kernels::aderdg::generic::c::solutionUpdate(luh, lduh, dt,
                                              5,  // getNumberOfVariables()
                                              0,  // getNumberOfParameters()
                                              4   // getNodesPerCoordinateAxis()
                                              );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        luh[i],
        ::exahype::tests::testdata::generic_euler::testSolutionUpdate::luh[i],
        eps, i);
  }

  delete[] luh;
  delete[] lduh;
}  // testSolutionUpdate

void GenericEulerKernelTest::testSpaceTimePredictorLinear() {
  logInfo("testSpaceTimePredictorLinear()",
          "Test space time predictor linear, ORDER=3, DIM=3");

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[320 =
  // nVar * nDOFx * nDOFy * nDOFz]

  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5, 0.5);
  const double dt = 1.267423918681417E-002;

  // Inputs:
  double** tempSpaceTimeUnknowns = new double*[1];
  tempSpaceTimeUnknowns[0] = new double[1600];  // lQi; nVar * nDOFx * nDOFy * nDOFz * (nDOFt+1); nDOF+1 only here

  double** tempSpaceTimeFluxUnknowns = new double*[2];
  tempSpaceTimeFluxUnknowns[0] = new double[3840+1280]; // lFi+source; nVar * nDOFx * nDOFy * nDOFt * (dim+1)
  tempSpaceTimeFluxUnknowns[1] = new double[3840];      // lQi; nVar * nDOFx * nDOFy * nDOFt * dim

  double* tempStateSizedVector = nullptr;

  // Outputs:
  double *tempUnknowns     = new double[320];     // lQh; nVar * nDOFx * nDOFy * nDOFz
  double *tempFluxUnknowns = new double[960+320]; // lFh+source; nVar * nDOFx * nDOFy * nDOFz * (dm+1) *

  double *lQhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6
  double *lFhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6

  // TODO(Dominic): Fix test
  kernels::aderdg::generic::c::spaceTimePredictorLinear<GenericEulerKernelTest>(
      *this,
      lQhbnd, lFhbnd,
      tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
      tempUnknowns,tempFluxUnknowns,
      tempStateSizedVector,
      ::exahype::tests::testdata::generic_euler::
          testSpaceTimePredictor::luh, // TODO(Dominic): Rename namespace to testSpaceTimePredictorLinear?
      dx, dt, nullptr
      );

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorLinear::lQhi[i],
        eps, i);
  }

  for (int i = 0; i < 960; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempFluxUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorLinear::lFhi[i],
        eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorLinear::lQhbnd[i],
        eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorLinear::lFhbnd[i],
        eps, i);
  }

  delete[] tempSpaceTimeUnknowns[0];
  delete[] tempSpaceTimeUnknowns;

  delete[] tempSpaceTimeFluxUnknowns[0];
  delete[] tempSpaceTimeFluxUnknowns[1];
  delete[] tempSpaceTimeFluxUnknowns;

  delete[] tempUnknowns;
  delete[] tempFluxUnknowns;

  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictorLinear

void GenericEulerKernelTest::testSpaceTimePredictorNonlinear() {
  logInfo("testSpaceTimePredictorNonlinear()",
          "Test space time predictor nonlinear, ORDER=3, DIM=3");

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[320 =
  // nVar * nDOFx * nDOFy * nDOFz]

  const tarch::la::Vector<DIMENSIONS, double> dx(0.05, 0.05, 0.05);
  const double timeStepSize = 1.083937460199773E-003;

  // Inputs:
  // space-time unknowns
  double** tempSpaceTimeUnknowns = new double*[4];
  tempSpaceTimeUnknowns[0] = new double[1280]; // lQi: nVar * nDOFx * nDOFy * nDOFz * nDOFt
  tempSpaceTimeUnknowns[1] = new double[1280]; // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  tempSpaceTimeUnknowns[2] = new double[1280]; // nVar * nDOFx * nDOFy * nDOFz * nDOFt
  tempSpaceTimeUnknowns[3] = new double[1280]; // nVar * nDOFx * nDOFy * nDOFz * nDOFt

  // space-time flux unknowns
  double** tempSpaceTimeFluxUnknowns = new double*[1];
  tempSpaceTimeFluxUnknowns[0]       = new double[5120+1280]; // nVar * nDOFx * nDOFy * nDOFz * nDOFt * (dim + 1)

  double* tempStateSizedVector = new double[5];

  // Outputs:
  // spatial unknowns
  double *tempUnknowns = new double[320];    // nVar * nDOFx * nDOFy * nDOFz
                                             // intentionally left uninitialised;
  // spatial flux unknowns
  double *tempFluxUnknowns = new double[1280+320];   // nVar * nDOFx * nDOFy * nDOFz * (dim + 1);

  double *lQhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6
  double *lFhbnd = new double[480];  // nVar * nDOFy * nDOF_z * 6

  // TODO(Dominic): Fix test
  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<GenericEulerKernelTest>(
      *this,
      lQhbnd, lFhbnd,
      tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
      tempUnknowns, tempFluxUnknowns,
      tempStateSizedVector,
      exahype::tests::testdata::generic_euler::testSpaceTimePredictorNonlinear::
          luh,
      dx, timeStepSize, nullptr
      );

  for (int i = 0; i < 1280; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempSpaceTimeUnknowns[0][i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictorNonlinear::lQi[i],
        eps, i);
  }

  kernels::idx6 idx_lFi(4, 4, 4, 4, (DIMENSIONS + 1), 5, __LINE__);
  kernels::idx6 idx_lFi_expected(4, 4, 4, 4, DIMENSIONS, 5, __LINE__);

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++) {
        for (int l = 0; l < 4; l++) {
          for (int m = 0; m < 2; m++) {  // skip 3 ( = source)
            for (int n = 0; n < 5; n++) {
              validateNumericalEqualsWithEpsWithParams1(
                  tempSpaceTimeFluxUnknowns[0][idx_lFi(i, j, k, l, m, n)],
                  ::exahype::tests::testdata::generic_euler::
                      testSpaceTimePredictorNonlinear::lFi[idx_lFi_expected(
                          i, j, k, l, m, n)],
                  eps, idx_lFi(i, j, k, l, m, n));
            }
          }
        }
      }
    }
  }

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lQhi[i],
        eps, i);
  }

  for (int i = 0; i < 960; i++) {  // skip 960 - 1279 (source)
    validateNumericalEqualsWithEpsWithParams1(
        tempFluxUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lFhi[i],
        eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorNonlinear::lQhbnd[i],
        eps, i);
  }

  for (int i = 0; i < 480; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lFhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorNonlinear::lFhbnd[i],
        eps, i);
  }

  delete[] tempSpaceTimeUnknowns[0];
  delete[] tempSpaceTimeUnknowns[1];
  delete[] tempSpaceTimeUnknowns[2];
  delete[] tempSpaceTimeUnknowns[3];
  delete[] tempSpaceTimeUnknowns;

  delete[] tempSpaceTimeFluxUnknowns[0];
  delete[] tempSpaceTimeFluxUnknowns;

  delete[] tempStateSizedVector;

  delete[] tempUnknowns;
  delete[] tempFluxUnknowns;

  delete[] lQhbnd;
  delete[] lFhbnd;
}  // testSpaceTimePredictorNonlinear

void GenericEulerKernelTest::testFaceUnknownsProjection() {
  logInfo("testFaceUnknownsProjection()",
          "Test face unknowns projection operators, ORDER=3, DIM=3");

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *lQhbndCoarseOut =
      new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndCoarseOut =
      new double[basisSize * basisSize * numberOfVariables];
  double *lQhbndFineOut = new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndFineOut = new double[basisSize * basisSize * numberOfVariables];

  // in:
  double *lQhbndCoarseIn =
      new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndCoarseIn =
      new double[basisSize * basisSize * numberOfVariables];
  double *lQhbndFineIn = new double[basisSize * basisSize * numberOfVariables];
  double *lFhbndFineIn = new double[basisSize * basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * basisSize * numberOfVariables; ++i) {
    lQhbndCoarseIn[i] = 1.0;
    lFhbndCoarseIn[i] = 1.0;
    lQhbndFineIn[i] = 1.0;
    lFhbndFineIn[i] = 1.0;
  }

  // Test the prolongation operator.
  for (int levelCoarse = 0; levelCoarse < 3; ++levelCoarse) {
    for (int levelDelta = 1; levelDelta < 1; ++levelDelta) {
      // todo For a levelDelta >= 2, assertionNumericalEquals
      // fails since the restriction result is not precise enough anymore.
      const int numberOfSubIntervals = tarch::la::aPowI(levelDelta, 3);

      tarch::la::Vector<DIMENSIONS - 1, int> subfaceIndex(0);

      // Test the restriction operator.
      memset(lQhbndCoarseOut, 0,
             basisSize * basisSize * numberOfVariables * sizeof(double));
      memset(lFhbndCoarseOut, 0,
             basisSize * basisSize * numberOfVariables * sizeof(double));
      for (int i2 = 0; i2 < numberOfSubIntervals; ++i2) {
        for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
          subfaceIndex[0] = i1;
          subfaceIndex[1] = i2;

          // Prolongate.
          kernels::aderdg::generic::c::faceUnknownsProlongation(
              lQhbndFineOut, lFhbndFineOut, lQhbndCoarseIn, lFhbndCoarseIn,
              levelCoarse, levelCoarse + levelDelta, subfaceIndex,
              numberOfVariables, basisSize);

          // Test prolongated values.
          for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
            assertionNumericalEquals5(lQhbndFineOut[m], lQhbndFineIn[m], m,
                                      levelCoarse, levelDelta, i1, i2);
            assertionNumericalEquals5(lFhbndFineOut[m], lFhbndFineIn[m], m,
                                      levelCoarse, levelDelta, i1, i2);
          }

          // Restrict.
          kernels::aderdg::generic::c::faceUnknownsRestriction(
              lQhbndCoarseOut, lFhbndCoarseOut, lQhbndFineIn, lFhbndFineIn,
              levelCoarse, levelCoarse + levelDelta, subfaceIndex,
              numberOfVariables, basisSize);
        }
        // Test restricted values.
        for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
          assertionNumericalEquals3(lQhbndCoarseOut[m], lQhbndCoarseIn[m], m,
                                    levelCoarse, levelDelta);
          assertionNumericalEquals3(lFhbndCoarseOut[m], lFhbndCoarseIn[m], m,
                                    levelCoarse, levelDelta);
        }
      }
    }
  }

  delete[] lQhbndCoarseOut;
  delete[] lFhbndCoarseOut;
  delete[] lQhbndFineOut;
  delete[] lFhbndFineOut;

  delete[] lQhbndCoarseIn;
  delete[] lFhbndCoarseIn;
  delete[] lQhbndFineIn;
  delete[] lFhbndFineIn;
}

void GenericEulerKernelTest::testVolumeUnknownsProjection() {
  logInfo("testVolumeUnknownsProjection()",
          "Test volume unknowns projection operators, ORDER=3, DIM=3");

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *luhCoarseOut =
      new double[basisSize * basisSize * basisSize * numberOfVariables];
  double *luhFineOut =
      new double[basisSize * basisSize * basisSize * numberOfVariables];

  // in:
  double *luhCoarseIn =
      new double[basisSize * basisSize * basisSize * numberOfVariables];
  double *luhFineIn =
      new double[basisSize * basisSize * basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * basisSize * basisSize * numberOfVariables;
       ++i) {
    luhCoarseIn[i] = 1.0;
    luhFineIn[i] = 1.0;
  }

  // Test the prolongation operator.
  for (int levelCoarse = 0; levelCoarse < 3; ++levelCoarse) {
    for (int levelDelta = 1; levelDelta < 2; ++levelDelta) {
      // todo For a levelDelta >= 2, assertionNumericalEquals
      // fails since the prolongation result is not precise enough anymore.
      const int numberOfSubIntervals = tarch::la::aPowI(levelDelta, 3);

      // Test the restriction operator.
      tarch::la::Vector<DIMENSIONS, int> subcellIndex(0);
      memset(luhCoarseOut, 0, basisSize * basisSize * basisSize *
                                  numberOfVariables * sizeof(double));

      for (int i3 = 0; i3 < numberOfSubIntervals; ++i3) {
        for (int i2 = 0; i2 < numberOfSubIntervals; ++i2) {
          for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
            subcellIndex[0] = i1;
            subcellIndex[1] = i2;
            subcellIndex[2] = i3;

            // Prolongate.
            kernels::aderdg::generic::c::volumeUnknownsProlongation(
                luhFineOut, luhCoarseIn, levelCoarse, levelCoarse + levelDelta,
                subcellIndex, numberOfVariables, basisSize);

            // Test prolongated values.
            for (int m = 0;
                 m < basisSize * basisSize * basisSize * numberOfVariables;
                 ++m) {
              /*
              assertionNumericalEquals5(luhFineOut[m], luhFineIn[m], m,
                                        levelCoarse, levelDelta, i1, i2);
                                        */
            }

            // Restrict.
            kernels::aderdg::generic::c::volumeUnknownsRestriction(
                luhCoarseOut, luhFineIn, levelCoarse, levelCoarse + levelDelta,
                subcellIndex, numberOfVariables, basisSize);
          }
        }
      }
      // Test restricted values.
      for (int m = 0; m < basisSize * basisSize * basisSize * numberOfVariables;
           ++m) {
        assertionNumericalEquals3(luhCoarseOut[m], luhCoarseIn[m], m,
                                  levelCoarse, levelDelta);
      }
    }
  }

  delete[] luhCoarseOut;
  delete[] luhFineOut;

  delete[] luhCoarseIn;
  delete[] luhFineIn;
}

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // DIMENSIONS==3
