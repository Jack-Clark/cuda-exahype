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

// @todo ExaHyPE Lizenz
#ifndef _EXAHYPE_TESTS_ELASTICITY_KERNEL_TEST_H_
#define _EXAHYPE_TESTS_ELASTICITY_KERNEL_TEST_H_

#include "peano/utils/Globals.h"
#include "tarch/logging/Log.h"
#include "tarch/tests/TestCase.h"

namespace exahype {
namespace tests {
namespace c {

static const int kNumberOfParameters = 3;
static const int kNumberOfVariables = 9 + kNumberOfParameters;
static const int kN = 4;
static const int kBasisSize = kN + 1;

class ElasticityKernelTest : public tarch::tests::TestCase {
 public:
  ElasticityKernelTest();
  virtual ~ElasticityKernelTest();

  void run() override;

 private:
  static tarch::logging::Log _log;

  //  void testPDEFluxes();
  void testSpaceTimePredictorLinear();
  //  void testSpaceTimePredictorNonlinear();
  void testVolumeIntegralLinear();
  //  void testVolumeIntegralNonlinear();
  void testRiemannSolverLinear();
  //  void testRiemannSolverNonlinear();
  void testSurfaceIntegralLinear();
  //  void testSurfaceIntegralNonlinear();
  //  void testSolutionUpdate();
  //  void testVolumeUnknownsProjection();
  //  void testFaceUnknownsProjection();
  //  void testEquidistantGridProjection();

 public:
  static int getNumberOfVariables() { return kNumberOfVariables; }
  static int getNumberOfParameters() { return kNumberOfParameters; }
  static int getNodesPerCoordinateAxis() { return kN; }

  static void flux(const double* const Q, double** F);

  static void source(const double* Q, double* S);

  static void eigenvalues(const double* const Q,
                              const int normalNonZeroIndex, double* lambda);

  static void ncp(const double* const Q, const double* const gradQ,
                      double* BgradQ);

  static void matrixb(const double* const Q, const int normalNonZero,
                          double* Bn);

  const double eps = 1.0e-9;  // for quick adaption of the test cases (say,
                              // switch to single precision)
};

}  // namespace c
}  // namespace tests
}  // namespace exahype

#endif  // _EXAHYPE_TESTS_ELASTICITY_KERNEL_TEST_H_
