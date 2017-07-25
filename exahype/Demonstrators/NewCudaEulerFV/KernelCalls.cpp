// This file is generated by the ExaHyPE toolkit.
// Please do not modify - it will be overwritten by the next
// ExaHyPE toolkit call.
// 
// ========================
//   www.exahype.eu
// ========================

#include <sstream>
#include <ostream>
#include "exahype/plotters/Plotter.h"
#include "exahype/profilers/ProfilerFactory.h"
#include "exahype/solvers/Solver.h"
#include "exahype/solvers/SolverCoupling.h"

#include "kernels/KernelCalls.h"

#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/GaussLobattoQuadrature.h"
#include "kernels/LimiterProjectionMatrices.h"
#include "kernels/DGMatrices.h"
#include "kernels/DGBasisFunctions.h"

#include "MyCudaEulerSolver.h"
#include "EulerWriter.h"



void kernels::initSolvers(exahype::Parser& parser, std::vector<std::string>& cmdlineargs) {
  {
  // Create and register solver
  exahype::solvers::RegisteredSolvers.push_back( new NewCudaEulerFV::MyCudaEulerSolver(parser.getMaximumMeshSize(0), parser.getTimeStepping(0)  , cmdlineargs));
  parser.checkSolverConsistency(0);

  
  }
  exahype::plotters::RegisteredPlotters.push_back( new exahype::plotters::Plotter(0,0,parser,new NewCudaEulerFV::EulerWriter(  *static_cast<NewCudaEulerFV::MyCudaEulerSolver*>(exahype::solvers::RegisteredSolvers[0])) ));


  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::initGaussLegendreNodesAndWeights(orders);
  kernels::initGaussLobattoNodesAndWeights(orders);
  kernels::initLimiterProjectionMatrices(orders);
  kernels::initDGMatrices(orders);
  kernels::initBasisFunctions(orders);
}


void kernels::finalise() {
  std::set<int> orders;
  for (const auto p : exahype::solvers::RegisteredSolvers) {
    orders.insert(p->getNodesPerCoordinateAxis()-1);
  }
  kernels::freeGaussLegendreNodesAndWeights(orders);
  kernels::freeGaussLobattoNodesAndWeights(orders);
  kernels::freeLimiterProjectionMatrices(orders);
  kernels::freeDGMatrices(orders);
  kernels::freeBasisFunctions(orders);

  for (auto solver : exahype::solvers::RegisteredSolvers) {
    delete solver;
  }
  exahype::solvers::RegisteredSolvers.clear();

  for (auto plotter : exahype::plotters::RegisteredPlotters) {
    delete plotter;
  }
  exahype::plotters::RegisteredPlotters.clear();
  for (auto coupling : exahype::solvers::RegisteredSolverCouplings) {
    delete coupling;
  }
  exahype::solvers::RegisteredSolverCouplings.clear();
}



void kernels::toString(std::ostream& ostream) {
/* Generated SolverRegistration code by the toolkit */

  ostream << "projectName: NewCudaEulerFV\n";
  ostream << "useOptimisedKernels: no\n";
  ostream << "Kernel[1].registration: FiniteVolumesSolver\n";
  ostream << "Kernel[1].type: ";
  exahype::solvers::RegisteredSolvers[1]->toString(ostream);
  ostream << "\n";
  ostream << "Kernel[1].hasConstants: null\n";
  ostream << "Kernel[1].patchSize: 10\n";
  ostream << "Kernel[0].Plotter[0]: NewCudaEulerFV::EulerWriter\n";
}


