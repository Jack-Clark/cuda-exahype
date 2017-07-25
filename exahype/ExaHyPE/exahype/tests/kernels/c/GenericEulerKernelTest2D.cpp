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

#include <cstring>

#include "../testdata/generic_euler_testdata.h"
#include "kernels/aderdg/generic/Kernels.h"


#if DIMENSIONS == 2

namespace exahype {
namespace tests {
namespace c {

// 2.5D Euler

void GenericEulerKernelTest::flux(const double *Q, double **F) {
  double *f = F[0];
  double *g = F[1];

  const double GAMMA = 1.4;

  const double irho = 1.0 / Q[0];
  const double p =
      (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

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
}

void GenericEulerKernelTest::source(const double *Q, double *S) {
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
  double p = (GAMMA - 1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

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
  // Arbitrary BS

  // Sven: I have no clue what these data shall do, but
  // ensure len(BgradQ)==len(Q).

  BgradQ[0] = 0;
  BgradQ[1] = 0;
  BgradQ[2] = 0;
  BgradQ[3] = 0;
  BgradQ[4] = 0;

  // Q[5]
  // gradQ[2][5]
  // BgradQ[2][5]
  if (!_setNcpAndMatrixBToZero) {
    BgradQ[0] = Q[0];
    BgradQ[1] = Q[3];
    BgradQ[2] = 3.0;
    BgradQ[3] = gradQ[0];
    BgradQ[4] = 0.7;
  }
}  // ncp

void GenericEulerKernelTest::matrixb(const double *const Q,
                                         const int normalNonZero, double *Bn) {
  std::memset(Bn, 0, 5 * 5 * sizeof(double));

  if (!_setNcpAndMatrixBToZero) {
    // 3D compressible Euler equations
    double *B1 = new double[5 * 5];
    double *B2 = new double[5 * 5];

    std::memset(B1, 0, 5 * 5 * sizeof(double));
    std::memset(B2, 0, 5 * 5 * sizeof(double));
    // Bn = B1 if normalNonZero == 0
    //      B2 if normalNonZero == 1
    std::memcpy(Bn, (normalNonZero == 0) ? B1 : B2, 5 * 5 * sizeof(double));

    delete[] B1;
    delete[] B2;
  }
}  // matrixb

void GenericEulerKernelTest::testPDEFluxes() {
  logInfo( "testPDEFluxes()", "Test PDE-related functions, DIM=2" );

  double Q[5] = {1., 0.1, 0.2, 0.3, 3.5};  // pressure = 1.39
  double f[5], g[5];
  double *F[2] = {f, g};

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

}  // testPDEFluxes

void GenericEulerKernelTest::testSolutionUpdate() {
  logInfo( "testSolutionUpdate()",  "Test solution update, ORDER=2, DIM=2" );

  // inputs:
  double *luh = new double[80];
  std::fill(luh, luh + 80, 0.0);
  for (int i = 0; i < 80; i += 5) {
    luh[i] = 1.0;
    luh[i + 4] = 2.5;
  }

  const double dt = 1.40831757919882352703e-03;
  // ::exahype::tests::testdata::generic_euler::GenericEulerKernelTest::lduh[80]

  kernels::aderdg::generic::c::solutionUpdate(
      luh, ::exahype::tests::testdata::generic_euler::testSolutionUpdate::lduh,
      dt,
      5,  // getNumberOfVariables(),
      0,  // getNumberOfParameters
      4   // getNodesPerCoordinateAxis()
      );

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        luh[i],
        ::exahype::tests::testdata::generic_euler::testSolutionUpdate::luh[i],
        eps, i);
  }

  delete[] luh;
}  // testSolutionUpdate

void GenericEulerKernelTest::testSurfaceIntegralLinear() {
  logInfo( "testSurfaceIntegralLinear()", "Test surface integral linear, ORDER=2, DIM=2" );

  {  // test 1
    // inputs:
    const double dx[2] = {0.1, 0.1};  // mesh spacing
    double lFhbnd[5 * 4 * 4] = {};    // nVar * nDofY * 4, zero initialized

    double *FLeft = &lFhbnd[0];
    double *FRight = &lFhbnd[20];
    double *FFront = &lFhbnd[40];
    double *FBack = &lFhbnd[60];

    for (int i = 0; i < 20; i += 5) {
      // in x orientation 1
      FLeft[i + 1] = 1.;
      FRight[i + 1] = 1.;
      // in y orientation 1
      FFront[i + 2] = 1.;
      FBack[i + 2] = 1.;
    }

    // input:
    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in
    double *lduh = new double[80];
    std::memcpy(
        lduh,
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in,
        80 * sizeof(double));

    // lFhbnd = [ FLeft | FRight | FFront | FBack ]
    kernels::aderdg::generic::c::surfaceIntegralLinear(
        lduh, lFhbnd, dx[0],
        5,  // getNumberOfVariables(),
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testSurfaceIntegralLinear::lduh_out_1[i],
          eps, i);
    }

    delete[] lduh;
  }  // test 1

  {  // test 2
    // inputs:
    const double dx[2] = {0.1, 0.1};  // mesh spacing
    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lFhbnd_in
    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in
    double *lduh = new double[80];
    std::memcpy(
        lduh,
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in,
        80 * sizeof(double));

    // lFhbnd = [ FLeft | FRight | FFront | FBack ]
    kernels::aderdg::generic::c::surfaceIntegralLinear(
        lduh, ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::
                  lFhbnd_in,
        dx[0],
        5,  // getNumberOfVariables(),
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testSurfaceIntegralLinear::lduh_out_2[i],
          eps, i);
    }

    delete[] lduh;
  }  // test 2

}  // testSurfaceIntegralLinear

void GenericEulerKernelTest::testSurfaceIntegralNonlinear() {
  logInfo( "testSurfaceIntegralNonlinear()", "Test surface integral nonlinear, ORDER=2, DIM=2" );

  {  // test 1
    // inputs:
    const double dx[2] = {0.1, 0.1};  // mesh spacing
    double lFhbnd[4 * 20] = {0.};

    double *FLeft = &lFhbnd[0];
    double *FRight = &lFhbnd[20];
    double *FFront = &lFhbnd[40];
    double *FBack = &lFhbnd[60];

    for (int i = 0; i < 20; i += 5) {
      // in x orientation 1
      FLeft[i + 1] = 1.;
      FRight[i + 1] = 1.;
      // in y orientation 1
      FFront[i + 2] = 1.;
      FBack[i + 2] = 1.;
    }

    // input:
    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in
    double *lduh = new double[80];
    std::memcpy(
        lduh,
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in,
        80 * sizeof(double));

    // lFhbnd = [ FLeft | FRight | FFront | FBack ]
    kernels::aderdg::generic::c::surfaceIntegralNonlinear(
        lduh, lFhbnd, dx[0],
        5,  // getNumberOfVariables(),
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testSurfaceIntegralNonlinear::lduh_1[i],
          eps, i);
    }

    delete[] lduh;
  }  // test 1

  {  // test 2
    // inputs:
    const double dx[2] = {0.1, 0.1};  // mesh spacing

    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in
    // ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lFhbnd_in

    double *lduh = new double[80];
    std::memcpy(
        lduh,
        ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::lduh_in,
        80 * sizeof(double));

    kernels::aderdg::generic::c::surfaceIntegralNonlinear(
        lduh, ::exahype::tests::testdata::generic_euler::testSurfaceIntegral::
                  lFhbnd_in,
        dx[0],
        5,  // getNumberOfVariables(),
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testSurfaceIntegralNonlinear::lduh_2[i],
          eps, i);
    }

    delete[] lduh;
  }  // test2
}  // testSurfaceIntegralNonlinear

void GenericEulerKernelTest::testRiemannSolverLinear() {
  // Rusanov
  logInfo( "testRiemannSolverLinear()", "Test Riemann Solver linear (Rusanov), ORDER=2, DIM=2" );

  {                                  // test normalNonZero = 0
                                     // output:
    double *FL = new double[5 * 4];  // nVar * nDOF(2)
    double *FR = new double[5 * 4];  // nVar * nDOF(2)
    // TODO: Improve test
    std::memset(FL, 0, 20 * sizeof(double));
    std::memset(FR, 0, 20 * sizeof(double));

    double  *tempFaceUnknowns      = new double[20]; // nDOF(1) * nVar
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

    // input
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR
    const int normalNonZero = 0;
    const double dt = 0.1;

    // TODO(Dominic): Fix test

    kernels::aderdg::generic::c::riemannSolverLinear<GenericEulerKernelTest>(
        *this,
        FL, FR,
        ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
        ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR,
         tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,
        dt,
        normalNonZero);

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverLinear::FL_1[i],
          eps, i);
    }

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverLinear::FR_1[i],
          eps, i);
    }

    delete[] FL;
    delete[] FR;
    delete[] tempStateSizedVectors[0];
    delete[] tempStateSizedVectors;
    delete[] tempStateSizedSquareMatrices[0];
    delete[] tempStateSizedSquareMatrices;
    delete[] tempFaceUnknowns;
  }  // end normalNonZero = 0

  {                                  // test normalNonZero = 1
                                     // output:
    double *FL = new double[5 * 4];  // nVar * nDOF(2)
    double *FR = new double[5 * 4];  // nVar * nDOF(2)
    // TODO: Improve test
    std::memset(FL, 0, 20 * sizeof(double));
    std::memset(FR, 0, 20 * sizeof(double));

    double  *tempFaceUnknowns      = new double[20]; // nDOF(1) * nVar
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

    // input:
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR
    const int normalNonZero = 1;
    const double dt = 0.1;

    // TODO(Dominic): Fix test
    kernels::aderdg::generic::c::riemannSolverLinear<GenericEulerKernelTest>(
        *this,
        FL, FR,
        ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL,
        ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR,
         tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,
        dt,
        normalNonZero
    );

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverLinear::FL_2[i],
          eps, i);
    }

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverLinear::FR_2[i],
          eps, i);
    }

    delete[] FL;
    delete[] FR;
    delete[] tempStateSizedVectors[0];
    delete[] tempStateSizedVectors;
    delete[] tempStateSizedSquareMatrices[0];
    delete[] tempStateSizedSquareMatrices;
    delete[] tempFaceUnknowns;
  }  // end normalNonZero = 1

}  // testRiemannSolverLinear

void GenericEulerKernelTest::testRiemannSolverNonlinear() {
  // Rusanov
  logInfo( "testRiemannSolverNonlinear()", "Test Riemann Solver nonlinear (Rusanov), ORDER=2, DIM=2" );

  {  // test 1
     // input:
    double QL[20] = {1., 0., 0., 0., 2.5, 1., 0., 0., 0., 2.5,
                     1., 0., 0., 0., 2.5, 1., 0., 0., 0., 2.5};
    double QR[20] = {1., 0., 0., 0., 2.5, 1., 0., 0., 0., 2.5,
                     1., 0., 0., 0., 2.5, 1., 0., 0., 0., 2.5};
    double  *tempFaceUnknownsArray        = nullptr;
    double **tempStateSizedVectors        = new double*[6];
    tempStateSizedVectors[0]              = new double[6*5];
    tempStateSizedVectors[1]              = tempStateSizedVectors[0]+5;
    tempStateSizedVectors[2]              = tempStateSizedVectors[1]+5;
    tempStateSizedVectors[3]              = tempStateSizedVectors[2]+5;
    tempStateSizedVectors[4]              = tempStateSizedVectors[3]+5;
    tempStateSizedVectors[5]              = tempStateSizedVectors[4]+5;
    double **tempStateSizedSquareMatrices = new double*[1];
    tempStateSizedSquareMatrices[0]       = new double[1*5*5];

    // inout:
    double *FL = new double[20];
    double *FR = new double[20];
    std::memset(FL, 0, 20 * sizeof(double));
    std::memset(FR, 0, 20 * sizeof(double));

    kernels::aderdg::generic::c::riemannSolverNonlinear<GenericEulerKernelTest>(
        *this,
        FL, FR, QL, QR,
        tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,
        0.0,  // dt
        0    // normalNonZero
        );

    // FL == FR, element by element
    for (int i = 0; i < 20; i++) {
      validateEquals(FL[i], FR[i]);
    }

    delete[] FL;
    delete[] FR;
    delete[] tempStateSizedSquareMatrices[0];
    delete[] tempStateSizedSquareMatrices;
    delete[] tempStateSizedVectors[0];
    delete[] tempStateSizedVectors;
  }  // test 1

  {  // test 2 nnz = 0
     // input:
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR

    // output:
    double *FL = new double[20];
    double *FR = new double[20];
    std::memcpy(FL, ::exahype::tests::testdata::generic_euler::
                        testRiemannSolverNonlinear::FL_1_in,
                20 * sizeof(double));
    std::memcpy(FR, ::exahype::tests::testdata::generic_euler::
                        testRiemannSolverNonlinear::FR_1_in,
                20 * sizeof(double));
    double  *tempFaceUnknownsArray        = nullptr;
    double **tempStateSizedVectors        = new double*[6];
    tempStateSizedVectors[0]              = new double[6*5];
    tempStateSizedVectors[1]              = tempStateSizedVectors[0]+5;
    tempStateSizedVectors[2]              = tempStateSizedVectors[1]+5;
    tempStateSizedVectors[3]              = tempStateSizedVectors[2]+5;
    tempStateSizedVectors[4]              = tempStateSizedVectors[3]+5;
    tempStateSizedVectors[5]              = tempStateSizedVectors[4]+5;
    double **tempStateSizedSquareMatrices = new double*[1];
    tempStateSizedSquareMatrices[0]       = new double[1*5*5];

    kernels::aderdg::generic::c::riemannSolverNonlinear<GenericEulerKernelTest>(
        *this,
        FL, FR, ::exahype::tests::testdata::generic_euler::
                    testRiemannSolverNonlinear::QL_1_in,
        ::exahype::tests::testdata::generic_euler::testRiemannSolverNonlinear::
            QR_1_in,
        tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,
        0.0,  // dt
        0     // normalNonZero
        );

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverNonlinear::F_1_out[i],
          eps, i);
    }

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverNonlinear::F_1_out[i],
          eps, i);
    }

    delete[] FL;
    delete[] FR;
    delete[] tempStateSizedSquareMatrices[0];
    delete[] tempStateSizedSquareMatrices;
    delete[] tempStateSizedVectors[0];
    delete[] tempStateSizedVectors;
  }  // test 2

  {  // test 3 nnz = 1
     // input:
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QL
    // ::exahype::tests::testdata::generic_euler::testRiemannSolver::QR

    // output:
    double *FL = new double[20];
    double *FR = new double[20];
    std::memcpy(FL, ::exahype::tests::testdata::generic_euler::
                        testRiemannSolverNonlinear::FL_2_in,
                20 * sizeof(double));
    std::memcpy(FR, ::exahype::tests::testdata::generic_euler::
                        testRiemannSolverNonlinear::FR_2_in,
                20 * sizeof(double));
    double  *tempFaceUnknownsArray        = nullptr;
    double **tempStateSizedVectors        = new double*[6];
    tempStateSizedVectors[0]              = new double[6*5];
    tempStateSizedVectors[1]              = tempStateSizedVectors[0]+5;
    tempStateSizedVectors[2]              = tempStateSizedVectors[1]+5;
    tempStateSizedVectors[3]              = tempStateSizedVectors[2]+5;
    tempStateSizedVectors[4]              = tempStateSizedVectors[3]+5;
    tempStateSizedVectors[5]              = tempStateSizedVectors[4]+5;
    double **tempStateSizedSquareMatrices = new double*[1];
    tempStateSizedSquareMatrices[0]       = new double[1*5*5];

    kernels::aderdg::generic::c::riemannSolverNonlinear<GenericEulerKernelTest>(
        *this,
        FL, FR, ::exahype::tests::testdata::generic_euler::
                    testRiemannSolverNonlinear::QL_2_in,
        ::exahype::tests::testdata::generic_euler::testRiemannSolverNonlinear::
            QR_2_in,
        tempFaceUnknownsArray,tempStateSizedVectors,tempStateSizedSquareMatrices,
        0.0,  // dt
        1    // normalNonZero
        );

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FL[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverNonlinear::F_2_out[i],
          eps, i);
    }

    for (int i = 0; i < 20; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          FR[i], ::exahype::tests::testdata::generic_euler::
                     testRiemannSolverNonlinear::F_2_out[i],
          eps, i);
    }

    delete[] FL;
    delete[] FR;
    delete[] tempStateSizedSquareMatrices[0];
    delete[] tempStateSizedSquareMatrices;
    delete[] tempStateSizedVectors[0];
    delete[] tempStateSizedVectors;
  }  // test3

}  // testRiemannSolverNonlinear

void GenericEulerKernelTest::testVolumeIntegralLinear() {
  logInfo( "testVolumeIntegralLinear()", "Test volume integral linear, ORDER=2, DIM=2" );

  {  // first test
    // output:
    double *lduh = new double[80];

    // input:
    double dx[2] = {3.70370370370370349811e-02,
                    3.70370370370370349811e-02};  // mesh spacing
    // ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lFhi[160]

    kernels::aderdg::generic::c::volumeIntegralLinear(
        lduh,
        ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lFhi,
        dx[0],
        5,  // getNumberOfVariables(),
        0,  // getNumberOfParameters
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testVolumeIntegralLinear::lduh_1[i],
          eps, i);
    }

    delete[] lduh;
  }  // scope limiter first test

  {  // second test, analogous to 3d seed

    // input:
    double dx[2] = {0.05, 0.05};     // mesh spacing
    double *lFhi = new double[160];  // nVar * dim * nDOFx * nDOFy
    std::fill(lFhi, lFhi + 160, 0.0);
    // lFhi = [ lFhi_x | lFhi_y ]
    double *lFhi_x = &lFhi[0];
    double *lFhi_y = &lFhi[80];

    // seed direction
    for (int i = 0; i < 80; i += 5) {
      lFhi_x[i + 1] = 1.;
      lFhi_y[i + 2] = 1.;
    }

    // output:
    double *lduh = new double[80];  // intentionally left uninitialised

    kernels::aderdg::generic::c::volumeIntegralLinear(lduh, lFhi, dx[0], 5, 0,
                                                      4);

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testVolumeIntegralLinear::lduh_2[i],
          eps, i);
    }

    delete[] lFhi;
    delete[] lduh;
  }  // scope second test
}  // testVolumeIntegralLinear

void GenericEulerKernelTest::testVolumeIntegralNonlinear() {
  logInfo( "testVolumeIntegralNonlinear()", "Test volume integral nonlinear, ORDER=2, DIM=2" );

  {  // first test

    // output:
    double *lduh = new double[80];

    // input:
    double dx[2] = {3.70370370370370349811e-02,
                    3.70370370370370349811e-02};  // mesh spacing
    // ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lFhi[240]

    kernels::aderdg::generic::c::volumeIntegralNonlinear(
        lduh,
        ::exahype::tests::testdata::generic_euler::testVolumeIntegral::lFhi,
        dx[0],
        5,  // getNumberOfVariables(),
        0,  // getNumberOfParameters()
        4   // getNodesPerCoordinateAxis()
        );

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testVolumeIntegralNonlinear::lduh_1[i],
          eps, i);
    }

    delete[] lduh;
  }  // scope limiter first test

  {  // second test, analogous to 3d seed

    // input:
    double dx[2] = {0.05, 0.05};     // mesh spacing
    double *lFhi = new double[240];  // nVar * (dim+1) * nDOFx * nDOFy
    std::fill(lFhi, lFhi + 240, 0.0);
    // lFhi = [ lFhi_x | lFhi_y | lShi]
    double *lFhi_x = &lFhi[0];
    double *lFhi_y = &lFhi[80];
    double *lShi = &lFhi[160];

    // seed direction
    for (int i = 0; i < 80; i += 5) {
      lFhi_x[i + 1] = 1.;
      lFhi_y[i + 2] = 1.;
    }

    std::fill(lShi, lShi + 80, 0.0);

    // output:
    double *lduh = new double[80];  // intentionally left uninitialised

    kernels::aderdg::generic::c::volumeIntegralNonlinear(lduh, lFhi, dx[0], 5,
                                                         0, 4);

    for (int i = 0; i < 80; i++) {
      validateNumericalEqualsWithEpsWithParams1(
          lduh[i], ::exahype::tests::testdata::generic_euler::
                       testVolumeIntegralNonlinear::lduh_2[i],
          eps, i);
    }

    delete[] lFhi;
    delete[] lduh;
  }  // scope limiter second test

}  // testVolumeIntegralNonlinear

void GenericEulerKernelTest::testSpaceTimePredictorLinear() {
  logInfo( "testSpaceTimePredictorLinear()", "Test space time predictor linear, ORDER=3, DIM=2" );

  // inputs:
  // exahype::tests::testdata::generic_euler::testSpaceTimePredictor::luh[80 =
  // nVar * nDOFx * nDOFy]

  const tarch::la::Vector<DIMENSIONS, double> dx(0.5, 0.5);
  const double dt = 1.267423918681417E-002;

  // Inputs:
  double** tempSpaceTimeUnknowns = new double*[1];
  tempSpaceTimeUnknowns[0] = new double[400];  // lQi; nVar * nDOFx * nDOFy * (nDOFt+1); nDOF+1 only here

  double** tempSpaceTimeFluxUnknowns = new double*[2];
  tempSpaceTimeFluxUnknowns[0] = new double[640+320];  // lFi; nDim * nVar * nDOFx * nDOFy * nDOFt
  tempSpaceTimeFluxUnknowns[1] = new double[640];  // gradQ; nDim * nVar * nDOFx * nDOFy * nDOFt

  double* tempStateSizedVector = nullptr;

  // Outputs:
  double *tempUnknowns     = new double[80];     // lQhi; nVar * nDOFx * nDOFz
  double *tempFluxUnknowns = new double[160+80]; // lFh,nVar * nDOFx * nDOFy * (dim+1)

  double *lQhbnd = new double[80];  // nVar * nDOFy * 4
  double *lFhbnd = new double[80];  // nVar * nDOFy * 4

  _setNcpAndMatrixBToZero = true;

  kernels::aderdg::generic::c::spaceTimePredictorLinear<GenericEulerKernelTest>(
      *this,
      lQhbnd, lFhbnd,
      tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
      tempUnknowns,tempFluxUnknowns,
      tempStateSizedVector,
      ::exahype::tests::testdata::generic_euler::
          testSpaceTimePredictorLinear::luh,
      dx, dt, nullptr);

  _setNcpAndMatrixBToZero = false;

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorLinear::lQhi[i],
        eps, i);
  }

  for (int i = 0; i < 160; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempFluxUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorLinear::lFhi[i],
        eps, i);
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                      testSpaceTimePredictorLinear::lQbnd[i],
        eps, i);
  }

  // TODO: The "fixed" kernel from the coding week doesn't compute lFbnd.

  /*
  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(lFbnd[i], 0.0, eps, i);
  }
  */
  delete[] lQhbnd;
  delete[] lFhbnd;

  delete[] tempSpaceTimeUnknowns[0];
  delete[] tempSpaceTimeUnknowns;

  delete[] tempSpaceTimeFluxUnknowns[0];
  delete[] tempSpaceTimeFluxUnknowns[1];
  delete[] tempSpaceTimeFluxUnknowns;

  delete[] tempUnknowns;
  delete[] tempFluxUnknowns;
}  // testSpaceTimePredictorLinear

void GenericEulerKernelTest::testSpaceTimePredictorNonlinear() {
  logInfo( "testSpaceTimePredictorNonlinear()", "Test space time predictor nonlinear, ORDER=3, DIM=2" );

  const tarch::la::Vector<DIMENSIONS, double> dx(5e-02, 5e-02);
  const double dt = 1.686854344081342E-003;

  // Inputs:
  double** tempSpaceTimeUnknowns = new double*[4];
  tempSpaceTimeUnknowns[0] = new double[320];  // lQi,nVar * nDOFx * nDOFy * nDOFt
  tempSpaceTimeUnknowns[1] = new double[320];  // nVar * nDOFx * nDOFy * nDOFt
  tempSpaceTimeUnknowns[2] = new double[320];  // nVar * nDOFx * nDOFy * nDOFt
  tempSpaceTimeUnknowns[3] = new double[320];  // nVar * nDOFx * nDOFy * nDOFt

  double** tempSpaceTimeFluxUnknowns = new double*[2];
  tempSpaceTimeFluxUnknowns[0] = new double[960+320]; // lFi+source;   nVar * nDOFx * nDOFy * nDOFt * (dim + 1)
  tempSpaceTimeFluxUnknowns[1] = new double[960];     // gradQ; nVar * nDOFx * nDOFy * nDOFt * dim

  double* tempStateSizedVector = new double[5];

  // Outputs:
  double *tempUnknowns     = new double[80];     // lQh; nVar * nDOFx * nDOFy
  double *tempFluxUnknowns = new double[240+80]; // lFh+source; nVar * nDOFx * nDOFy * (dim + 1)

  double *lQhbnd = new double[80];  // nVar * nDOFy * 4
  double *lFhbnd = new double[80];  // nVar * nDOFy * 4

  _setNcpAndMatrixBToZero = true;

  kernels::aderdg::generic::c::spaceTimePredictorNonlinear<GenericEulerKernelTest>(
      *this,
      lQhbnd, lFhbnd,
      tempSpaceTimeUnknowns,tempSpaceTimeFluxUnknowns,
      tempUnknowns,tempFluxUnknowns,
      tempStateSizedVector,
      ::exahype::tests::testdata::generic_euler::
          testSpaceTimePredictorNonlinear::luh,
      dx, dt, nullptr);

  _setNcpAndMatrixBToZero = false;

  for (int i = 0; i < 320; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempSpaceTimeUnknowns[0][i], ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictorNonlinear::lQi[i],
        eps, i);
  }

  kernels::idx5 idx_lFi(4, 4, 4, (DIMENSIONS + 1), 5, __LINE__);
  kernels::idx5 idx_lFi_expected(4, 4, 4, DIMENSIONS, 5, __LINE__);

  for (int i = 0; i < 4; i++) {
    for (int k = 0; k < 4; k++) {
      for (int l = 0; l < 4; l++) {
        for (int m = 0; m < 2; m++) {  // skip 2 ( = source)
          for (int n = 0; n < 5; n++) {
            validateNumericalEqualsWithEpsWithParams1(
                tempSpaceTimeFluxUnknowns[0][idx_lFi(i, k, l, m, n)],
                ::exahype::tests::testdata::generic_euler::
                    testSpaceTimePredictorNonlinear::lFi[idx_lFi_expected(
                        i, k, l, m, n)],
                eps, idx_lFi(i, k, l, m, n));
          }
        }
      }
    }
  }

  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        tempUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lQhi[i],
        eps, i);
  }

  for (int i = 0; i < 160; i++) {  // skip 160 - 239 (source)
    validateNumericalEqualsWithEpsWithParams1(
        tempFluxUnknowns[i], ::exahype::tests::testdata::generic_euler::
                     testSpaceTimePredictorNonlinear::lFhi[i],
        eps, i);
  }

  // lQhbnd
  for (int i = 0; i < 80; i++) {
    validateNumericalEqualsWithEpsWithParams1(
        lQhbnd[i], ::exahype::tests::testdata::generic_euler::
                       testSpaceTimePredictorNonlinear::lQhbnd[i],
        eps, i);
  }

  // lFhbnd
  for (int i = 0; i < 80; i++) {
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
  logInfo( "testFaceUnknownsProjection()", "Test face unknowns projection operators, ORDER=2, DIM=2" );

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *lQhbndCoarseOut = new double[basisSize * numberOfVariables];
  double *lFhbndCoarseOut = new double[basisSize * numberOfVariables];
  double *lQhbndFineOut = new double[basisSize * numberOfVariables];
  double *lFhbndFineOut = new double[basisSize * numberOfVariables];

  // in:
  double *lQhbndCoarseIn = new double[basisSize * numberOfVariables];
  double *lFhbndCoarseIn = new double[basisSize * numberOfVariables];
  double *lQhbndFineIn = new double[basisSize * numberOfVariables];
  double *lFhbndFineIn = new double[basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * numberOfVariables; ++i) {
    lQhbndCoarseIn[i] = 1.0;
    lFhbndCoarseIn[i] = 1.0;
    lQhbndFineIn[i] = 1.0;
    lFhbndFineIn[i] = 1.0;
  }

  // Test the prolongation operator.
  for (int levelCoarse = 0; levelCoarse < 3; ++levelCoarse) {
    for (int levelDelta = 1; levelDelta < 3; ++levelDelta) {
      // todo For a levelDelta >= 4, assertionNumericalEquals
      // fails since the restriction result is not precise enough anymore.
      const int numberOfSubIntervals = tarch::la::aPowI(levelDelta, 3);

      // Test the restriction operator.
      memset(lQhbndCoarseOut, 0,
             sizeof(double) * numberOfVariables * basisSize);
      memset(lFhbndCoarseOut, 0,
             sizeof(double) * numberOfVariables * basisSize);
      for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
        // Prolongate.
        tarch::la::Vector<DIMENSIONS - 1, int> subfaceIndex(i1);
        kernels::aderdg::generic::c::faceUnknownsProlongation(
            lQhbndFineOut, lFhbndFineOut, lQhbndCoarseIn, lFhbndCoarseIn,
            levelCoarse, levelCoarse + levelDelta, subfaceIndex,
            numberOfVariables, basisSize);

        // Test prolongated values.
        for (int m = 0; m < basisSize * numberOfVariables; ++m) {
          assertionNumericalEquals(lQhbndFineOut[m], lQhbndFineIn[m]);
          assertionNumericalEquals(lFhbndFineOut[m], lFhbndFineIn[m]);
        }

        // Restrict.
        kernels::aderdg::generic::c::faceUnknownsRestriction(
            lQhbndCoarseOut, lFhbndCoarseOut, lQhbndFineIn, lFhbndFineIn,
            levelCoarse, levelCoarse + levelDelta, subfaceIndex,
            numberOfVariables, basisSize);
      }
      // Test restricted values.
      for (int m = 0; m < basisSize * numberOfVariables; ++m) {
        assertionNumericalEquals(lQhbndCoarseOut[m], lQhbndCoarseIn[m]);
        assertionNumericalEquals(lFhbndCoarseOut[m], lFhbndCoarseIn[m]);
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
  logInfo( "testVolumeUnknownsProjection()", "Test volume unknowns projection operators, ORDER=2, DIM=2" );

  const int numberOfVariables = 1;
  const int basisSize = 4;

  // in/out
  double *luhCoarseOut = new double[basisSize * basisSize * numberOfVariables];
  double *luhFineOut = new double[basisSize * basisSize * numberOfVariables];

  // in:
  double *luhCoarseIn = new double[basisSize * basisSize * numberOfVariables];
  double *luhFineIn = new double[basisSize * basisSize * numberOfVariables];

  // Initialise to constant value.
  for (int i = 0; i < basisSize * basisSize * numberOfVariables; ++i) {
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
      memset(luhCoarseOut, 0,
             basisSize * basisSize * numberOfVariables * sizeof(double));

      for (int i2 = 0; i2 < numberOfSubIntervals; ++i2) {
        for (int i1 = 0; i1 < numberOfSubIntervals; ++i1) {
          subcellIndex[0] = i1;
          subcellIndex[1] = i2;

          // Prolongate.
          kernels::aderdg::generic::c::volumeUnknownsProlongation(
              luhFineOut, luhCoarseIn, levelCoarse, levelCoarse + levelDelta,
              subcellIndex, numberOfVariables, basisSize);

          // Test prolongated values.
          for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
            assertionNumericalEquals5(luhFineOut[m], luhFineIn[m], m,
                                      levelCoarse, levelDelta, i1, i2);
          }

          // Restrict.
          kernels::aderdg::generic::c::volumeUnknownsRestriction(
              luhCoarseOut, luhFineIn, levelCoarse, levelCoarse + levelDelta,
              subcellIndex, numberOfVariables, basisSize);
        }
      }
      // Test restricted values.
      for (int m = 0; m < basisSize * basisSize * numberOfVariables; ++m) {
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

#endif  // DIMENSIONS==2
