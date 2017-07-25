#include "MyCudaEulerSolver.h"
#include "MyCudaEulerSolver_Variables.h"
#include "Logo.h"

/*extern "C" {
  void cudaEigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda);
  void cudaFlux(const double* const Q, double** F);
  void cudaFlux2D(const double* const Q, double** F);
}*/

void NewCudaEulerFV::MyCudaEulerSolver::init(std::vector<std::string>& cmdlineargs) {
  // @todo Please implement/augment if required
}

bool NewCudaEulerFV::MyCudaEulerSolver::hasToAdjustSolution(const tarch::la::Vector<DIMENSIONS, double>& center, const tarch::la::Vector<DIMENSIONS, double>& dx, const double t, const double dt) {
  // @todo Please implement/augment if required
  return tarch::la::equals(t,0.0);
}

void NewCudaEulerFV::MyCudaEulerSolver::adjustedSolutionValues(const double* const x,const double w,const double t,const double dt, double* Q) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  tarch::la::Vector<DIMENSIONS,double> myX( x[0] - 0.06, 1.0-x[1] - 0.25 ); // translate
  myX *= static_cast<double>(Image.width);
  tarch::la::Vector<DIMENSIONS,int>    myIntX( 1.2*myX(0) , 1.2*myX(1) );  // scale

  double Energy = 0.1;

  if (
    myIntX(0) > 0 && myIntX(0) < static_cast<int>(Image.width)
    &&
    myIntX(1) > 0 && myIntX(1) < static_cast<int>(Image.height)
  ) {
    Energy += 1.0-Image.pixel_data[myIntX(1)*Image.width+myIntX(0)];
  }

  Q[0] = 1.0;
  Q[1] = 0.0;
  Q[2] = 0.0;
  Q[3] = 0.0;
  Q[4] = Energy;

}

exahype::solvers::Solver::RefinementControl NewCudaEulerFV::MyCudaEulerSolver::refinementCriterion(const double* luh, const tarch::la::Vector<DIMENSIONS, double>& center,const tarch::la::Vector<DIMENSIONS, double>& dx, double t,const int level) {
  // @todo Please implement/augment if required
  return exahype::solvers::Solver::RefinementControl::Keep;
}


void NewCudaEulerFV::MyCudaEulerSolver::eigenvalues(const double* const Q, const int normalNonZeroIndex, double* lambda) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);   

  const double u_n = Q[normalNonZeroIndex + 1] * irho;
  const double c = sqrt(GAMMA * p * irho);

  lambda[0] = u_n - c;
  lambda[1] = u_n;
  lambda[2] = u_n;
  lambda[3] = u_n;
  lambda[4] = u_n + c;
}

void NewCudaEulerFV::MyCudaEulerSolver::flux(const double* const Q, double** F) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  const double GAMMA = 1.4;

  const double irho = 1.0/Q[0];
  const double p = (GAMMA-1) * (Q[4] - 0.5 * (Q[1] * Q[1] + Q[2] * Q[2]) * irho);

  double* f = F[0];
  double* g = F[1];

  f[0] = Q[1]; 
  f[1] = irho * Q[1] * Q[1] + p;
  f[2] = irho * Q[1] * Q[2];
  f[3] = irho * Q[1] * Q[3];
  f[4] = irho * Q[1] * (Q[4] + p);

  g[0] = Q[2];
  g[1] = irho * Q[1] * Q[2];
  g[2] = irho * Q[2] * Q[2] + p;
  g[3] = irho * Q[2] * Q[3];
  g[4] = irho * Q[2] * (Q[4] + p);
}


void NewCudaEulerFV::MyCudaEulerSolver::source(const double* const Q, double* S) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)
  
  S[0] = 0.0;
  S[1] = 0.0;
  S[2] = 0.0;
  S[3] = 0.0;
  S[4] = 0.0;
}

void NewCudaEulerFV::MyCudaEulerSolver::boundaryValues(
    const double* const x,
    const double t,const double dt,
    const int faceIndex,
    const int normalNonZero,
    const double* const stateInside,
    double* stateOutside) {
  // Dimensions             = 2
  // Number of variables    = 5 (#unknowns + #parameters)

  stateOutside[0] = stateInside[0];
  stateOutside[1] = stateInside[1];
  stateOutside[2] = stateInside[2];
  stateOutside[3] = stateInside[3];
  stateOutside[4] = stateInside[4];
}