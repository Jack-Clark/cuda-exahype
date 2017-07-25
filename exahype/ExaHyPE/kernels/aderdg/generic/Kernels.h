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

#ifndef _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_

#include "string.h"

#include "peano/utils/Globals.h"
#include "tarch/la/Vector.h"

#include "kernels/GaussLegendreQuadrature.h"

#include "kernels/DGMatrices.h"

#include "tarch/la/Scalar.h"
#include "tarch/la/ScalarOperations.h"

#include "tarch/logging/Log.h"

#define MbasisSize 4
#define Mvar 5
#define Mdim 3
#define f2p5(var, dim, i, j, k)                                        \
  (var + Mvar * dim + Mvar * Mdim * i + Mvar * Mdim * MbasisSize * j + \
   Mvar * Mdim * MbasisSize * MbasisSize * k)
#define p2f5(var, dim, i, j, k)                                   \
  (dim * MbasisSize * MbasisSize * MbasisSize * Mvar + Mvar * i + \
   Mvar * MbasisSize * j + Mvar * MbasisSize * MbasisSize * k + var)

#define Mface 6
#define f2p4(var, face, a, b) \
  (var + Mvar * face + Mvar * Mface * a + Mvar * Mface * MbasisSize * b)
#define p2f4(var, face, a, b)                                                 \
  (face * MbasisSize * MbasisSize * Mvar + Mvar * a + Mvar * MbasisSize * b + \
   var)

/** Computes a 1-d node index.
  * The brackets around \p ix allow to write
  * idx1(ix+ox,...), where ox is some offset.
  *
  * TODO(Dominic): Get rid of this; only used in c/2d/boundaryConditions.cpph
  */
#define nidx1(ix) numberOfVariables*(ix)

/** Computes a 2-d node index.
  * The brackets around \p ix and \p iy allow to write
  * idx1(ix+ox, iy+oy,...), where ox and oy are some offsets.
  *
  * TODO(Dominic): Get rid of this; only used in c/3d/boundaryConditions.cpph
  */
#define nidx2(ix, iy) numberOfVariables*(basisSize * (iy) + ix)

// todo Dominic Etienne Charrier
// Possibly redundant definition of face indices
// see exahype/solvers/Solver.h
// On the other hand, the kernels should be
// more or less independent of ExaHyPE/exahype.
#define EXAHYPE_FACE_LEFT 0
#define EXAHYPE_FACE_RIGHT 1
#define EXAHYPE_FACE_FRONT 2
#define EXAHYPE_FACE_BACK 3
#define EXAHYPE_FACE_BOTTOM 4
#define EXAHYPE_FACE_TOP 5

namespace kernels {
namespace aderdg {
namespace generic {
namespace c {

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
void spaceTimePredictorNonlinear(
    SolverType& solver,
    double*  lQhbnd, double* lFhbnd,
    double** tempSpaceTimeUnknowns,
    double** tempSpaceTimeFluxUnknowns,
    double*  tempUnknowns,
    double*  tempFluxUnknowns,
    double*  tempStateSizedVector,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double dt,
    double* tempPointForceSources);

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int numberOfParameters,
                    const int basisSize);

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                          const tarch::la::Vector<DIMENSIONS, double>& dx,
                          const int numberOfVariables,
                          const int numberOfParameters, const int basisSize);

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables,
                             const int numberOfParameters, const int basisSize);

// todo 10/02/16: Dominic
// Keep only one surfaceIntegral.
void surfaceIntegralNonlinear(double* lduh, const double* const lFbnd,
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const int numberOfVariables, const int basisSize);

void surfaceIntegralLinear(double* lduh, const double* const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double>& dx,
                           const int numberOfVariables, const int basisSize);

/*void surfaceIntegral2(
    double* lduh,
    const double* const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>&  dx,
    const int numberOfVariables,
    const int basisSize
);*/

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
void solutionAdjustment(SolverType& solver, double* luh,
                        const tarch::la::Vector<DIMENSIONS, double>& center,
                        const tarch::la::Vector<DIMENSIONS, double>& dx,
                        const double t, const double dt);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments
// template argument functions and non-template argument function.
template <typename SolverType>
void riemannSolverNonlinear(
    SolverType& solver,
    double* FL, double* FR, const double* const QL,
    const double* const QR,
    double*  tempFaceUnknowns,
    double** tempStateSizedVectors,
    double** tempStateSizedSquareMatrices,
    const double dt,
    const int normalNonZero);

template <typename SolverType>
void boundaryConditions(
    SolverType& solver,
    double* fluxOut, double* stateOut,
    const double* const fluxIn, const double* const stateIn, // TODO(Dominic): Inconsistent order of arguments w.r.t to riemannSolver
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize,
    const double t, const double dt, const int faceIndex,
    const int normalNonZero);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
double stableTimeStepSize(SolverType& solver, const double* const luh,
                          double* tempEigenvalues,
                          const tarch::la::Vector<DIMENSIONS, double>& dx);

void faceUnknownsProlongation(
    double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
    const double* lFhbndCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void faceUnknownsRestriction(
    double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
    const double* lFhbndFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsProlongation(
    double* luhFine, const double* luhCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsRestriction(
    double* luhCoarse, const double* luhFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);
    
//TODO KD    
template <typename SolverType>
void dummyK_Kernel(
    SolverType& solver,
    const double t,
    const double dt,
    const tarch::la::Vector<DIMENSIONS, double>& center,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const int numberOfVariables, 
    const int numberOfParameters, 
    const int basisSize,
    double* tempPointForceSources //memory space for forceVector
    );
}  // namespace c
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels

#if DIMENSIONS == 2
#include "kernels/aderdg/generic/c/2d/boundaryConditions.cpph"
#include "kernels/aderdg/generic/c/2d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/c/2d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/c/2d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/c/2d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/c/2d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/c/2d/stableTimeStepSize.cpph"
#include "kernels/aderdg/generic/c/2d/dummyK_Kernel.cpph"
#elif DIMENSIONS == 3
#include "kernels/aderdg/generic/c/3d/boundaryConditions.cpph"
#include "kernels/aderdg/generic/c/3d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/c/3d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/c/3d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/c/3d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/c/3d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/c/3d/stableTimeStepSize.cpph"
#include "kernels/aderdg/generic/c/3d/dummyK_Kernel.cpph"
#endif

// Todo: Recasting the code from function templates to class templates
//       did not yet consider the Fortran kernels and probably never will,
//       as they are different anyway.

namespace kernels {
namespace aderdg {
namespace generic {
namespace fortran {

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
void spaceTimePredictorNonlinear(
    SolverType& solver,
    double*  lQhbnd, double* lFhbnd,
    double** tempSpaceTimeUnknowns,
    double** tempSpaceTimeFluxUnknowns,
    double*  tempUnknowns,
    double*  tempFluxUnknowns,
    double*  tempStateSizedVector,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double dt,
    double* tempPointForceSources);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
void spaceTimePredictorLinear(
    SolverType& solver,
    double*  lQhbnd, double* lFhbnd,
    double** tempSpaceTimeUnknowns,
    double** tempSpaceTimeFluxUnknowns,
    double*  tempUnknowns,
    double*  tempFluxUnknowns,
    double*  tempStateSizedVector,
    const double* const luh,
    const tarch::la::Vector<DIMENSIONS, double>& dx,
    const double dt,
    double* tempPointForceSources);

/**
 * (At the moment, we always evaluate the time averaged space-time
 * predictor unknowns.)
 * todo docu
 */
void predictor(double* lQhi, double* lFhi, const double* const lQi,
               const double* const lFi, const double predictorTimeStepSize,
               const int numberOfVariables, const int basisSize);

/**
 * @todo Dominic Etienne Charrier
 * This is just a "parent" function that
 * invokes the function going by the same
 * name 2*dim times.
 */
void extrapolatedPredictor(double* lQhbnd, double* lFhbnd,
                           const double* const lQhi, const double* const lFhi,
                           const double predictorTimeStepSize,
                           const int numberOfVariables, const int basisSize);

// todo Dominic Etienne Charrier:
// The DIMENSIONS depending mesh size vector enables overloading at the moment.
// If we replace it by scalar mesh size, we have to add a template argument "int
// dim".

void solutionUpdate(double* luh, const double* const lduh, const double dt,
                    const int numberOfVariables, const int numberOfParameters, 
                    const int basisSize);

void volumeIntegralNonlinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables, const int numberOfParameters, 
                             const int basisSize);

void volumeIntegralLinear(double* lduh, const double* const lFhi,
                             const tarch::la::Vector<DIMENSIONS, double>& dx,
                             const int numberOfVariables, const int numberOfParameters, 
                             const int basisSize);

// todo 10/02/16: Dominic
// Keep only one surfaceIntegral.
void surfaceIntegralNonlinear(double* lduh, const double* const lFbnd,
                              const tarch::la::Vector<DIMENSIONS, double>& dx,
                              const int numberOfVariables, const int basisSize);

void surfaceIntegralLinear(double* lduh, const double* const lFbnd,
                           const tarch::la::Vector<DIMENSIONS, double>& dx,
                           const int numberOfVariables, const int basisSize);

/*void surfaceIntegral2(
    double* lduh,
    const double* const lFhbnd,
    const tarch::la::Vector<DIMENSIONS,double>&  dx,
    const int numberOfVariables,
    const int basisSize
);*/

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
void solutionAdjustment(SolverType& solver, double* luh,
                        const tarch::la::Vector<DIMENSIONS, double>& center,
                        const tarch::la::Vector<DIMENSIONS, double>& dx,
                        const double t, const double dt);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments
// template argument functions and non-template argument function.
template <typename SolverType>
void riemannSolverNonlinear(
    SolverType& solver,
    double* FL, double* FR, const double* const QL,
    const double* const QR,
    double*  tempFaceUnknowns,
    double** tempStateSizedVectors,
    double** tempStateSizedSquareMatrices,
    const double dt,
    const int normalNonZero);

template <typename SolverType>
void riemannSolverLinear(
    SolverType& solver,
    double* FL, double* FR,
    const double* const QL, const double* const QR,
    double*  tempFaceUnknowns,
    double** tempStateSizedVectors,
    double** tempStateSizedSquareMatrices,
    const double dt,
    const int normalNonZero);

// @todo Dominic Etienne Charrier
// Inconsistent ordering of inout and in arguments for
// template argument functions and non-template argument function.
template <typename SolverType>
double stableTimeStepSize(SolverType& solver, const double* const luh,
                          double* tempEigenvalues,
                          const tarch::la::Vector<DIMENSIONS, double>& dx);

void faceUnknownsProlongation(
    double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
    const double* lFhbndCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void faceUnknownsRestriction(
    double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
    const double* lFhbndFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsProlongation(
    double* luhFine, const double* luhCoarse, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

void volumeUnknownsRestriction(
    double* luhCoarse, const double* luhFine, const int coarseGridLevel,
    const int fineGridLevel,
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int numberOfVariables, const int basisSize);

}  // namespace fortran
}  // namespace generic
}  // namespace aderdg
}  // namespace kernels
#if DIMENSIONS == 3
#include "kernels/aderdg/generic/fortran/3d/riemannSolverLinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/riemannSolverNonlinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/solutionAdjustment.cpph"
#include "kernels/aderdg/generic/fortran/3d/spaceTimePredictorLinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/spaceTimePredictorNonlinear.cpph"
#include "kernels/aderdg/generic/fortran/3d/stableTimeStepSize.cpph"
// #elif DIMENSIONS == 2
// //@todo
// #include "kernels/aderdg/generic/fortran/2d/solutionAdjustment.cpph"
// #include "kernels/aderdg/generic/fortran/2d/stableTimeStepSize.cpph"
// #include "kernels/aderdg/generic/fortran/2d/spaceTimePredictorNonlinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/spaceTimePredictorLinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/riemannSolverNonlinear.cpph"
// #include "kernels/aderdg/generic/fortran/2d/riemannSolverLinear.cpph"
#endif

#endif /* _EXAHYPE_KERNELS_ADERDG_GENERIC_PDEFLUXES_H_ */
