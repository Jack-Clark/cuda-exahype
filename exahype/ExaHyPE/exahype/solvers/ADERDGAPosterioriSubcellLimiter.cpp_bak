/*
 * ADERDGAPosterioriSubcellLimiter.cpp
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#include "ADERDGAPosterioriSubcellLimiter.h"

#include "tarch/Assertions.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "kernels/limiter/generic/Limiter.h"

tarch::logging::Log exahype::solvers::ADERDGAPosterioriSubcellLimiter::_log("exahype::solvers::ADERDGAPosterioriSubcellLimiter");

exahype::solvers::ADERDGAPosterioriSubcellLimiter::ADERDGAPosterioriSubcellLimiter(
    int aderdgSolverNumber,int finiteVolumesSolverNumber)
      : CellWiseCoupling(0.0,0.0),
        _aderdgSolverNumber(aderdgSolverNumber),
        _finiteVolumesSolverNumber(finiteVolumesSolverNumber)
{
  assertion2(_aderdgSolverNumber >= 0 && _aderdgSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _aderdgSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion2(_finiteVolumesSolverNumber >= 0 && _finiteVolumesSolverNumber<static_cast<int>(exahype::solvers::RegisteredSolvers.size()),
      _finiteVolumesSolverNumber,exahype::solvers::RegisteredSolvers.size());
  assertion1(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()==exahype::solvers::Solver::Type::ADER_DG,
      exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]->getType()));
  assertion1(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()==exahype::solvers::Solver::Type::FiniteVolumes,
        exahype::solvers::Solver::toString(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]->getType()));
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::coupleFirstTime(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::ADERDGSolver* aderdgSolver =
      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
  exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
      static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
  int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
  int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);

  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
    auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[aderdgElement];
    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[finiteVolumesElement];

    // Set the initial conditions of the finite volumes solver.
    finiteVolumesSolver->synchroniseTimeStepping(cellDescriptionsIndex,finiteVolumesElement);

    if (finiteVolumesSolver->hasToAdjustSolution(
        finiteVolumesCellDescription.getOffset()+0.5*finiteVolumesCellDescription.getSize(),
        finiteVolumesCellDescription.getSize(),finiteVolumesCellDescription.getTimeStamp())) {
      finiteVolumesSolver->setInitialConditions(
          cellDescriptionsIndex,finiteVolumesElement,
          fineGridVertices,fineGridVerticesEnumerator);

      // Set the initial conditions of the ADER-DG solver by projecting
      // the solution values of the FV solver onto the DG solution space.
      aderdgSolver->synchroniseTimeStepping(cellDescriptionsIndex,aderdgElement);
      double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
          getData(finiteVolumesCellDescription.getSolution()).data();
      double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
          getData(aderdgCellDescription.getSolution()).data();
      kernels::limiter::generic::c::projectOnADERDGSpace(
          finiteVolumesSolution,
          finiteVolumesSolver->getNumberOfVariables(),
          finiteVolumesSolver->getNodesPerCoordinateAxis(),
          aderdgSolver->getNodesPerCoordinateAxis(),
          aderdgSolution);

      // This writes the min max values to face 0;
      double* solutionMin = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMin()).data();
      double* solutionMax = DataHeap::getInstance().getData(aderdgCellDescription.getSolutionMax()).data();
      std::fill(solutionMin,solutionMin+aderdgSolver->getNumberOfVariables(),std::numeric_limits<double>::max());
      std::fill(solutionMax,solutionMax+aderdgSolver->getNumberOfVariables(),-std::numeric_limits<double>::max()); // !!! Mind the "-"
      // kernels::limiter::generic::c::findCellLocalLimMinAndMax( //TODO Dominic: deprecated method
          // finiteVolumesSolution,
          // aderdgSolver->getNumberOfVariables(),
          // aderdgSolver->getNodesPerCoordinateAxis(),
          // solutionMin,
          // solutionMax);
      // This writes the face 0 min max values to the other faces
      for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
        std::copy(solutionMin, solutionMin+aderdgSolver->getNumberOfVariables(),
            solutionMin+i*aderdgSolver->getNumberOfVariables());
        std::copy(solutionMax, solutionMax+aderdgSolver->getNumberOfVariables(),
            solutionMax+i*aderdgSolver->getNumberOfVariables());
      }

      // Check if initial solution is troubled
      bool initialSolutionIsTroubled =
      false; //TODO Dominic: isTroubled now anticipates the new DG solution
//          kernels::limiter::generic::c::isTroubledCell(
//              aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
//              solutionMin,solutionMax);

      if (initialSolutionIsTroubled) {
        aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled); // implicit conversion.
      }
    }
  }
}

bool
oneNeighbourIsTroubled(const tarch::la::Vector<DIMENSIONS_TIMES_TWO,
    exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus>& limiterStatus) {
  for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
    switch (limiterStatus[i]) {
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell:
      return true;
      break;
    default:
      break;
    }
  }

  return false;
}

void evolveFiniteVolumesSolution(
    const int cellDescriptionsIndex,
    const int finiteVolumesElement,
    const int aderdgElement,
    exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver,
    exahype::solvers::ADERDGSolver* aderdgSolver,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
      cellDescriptionsIndex)[finiteVolumesElement];
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
      cellDescriptionsIndex)[aderdgElement];

  finiteVolumesCellDescription.setTimeStamp(aderdgCellDescription.getCorrectorTimeStamp());
  finiteVolumesCellDescription.setTimeStepSize(aderdgCellDescription.getCorrectorTimeStepSize());
  finiteVolumesSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

  double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
      getData(finiteVolumesCellDescription.getSolution()).data();
  double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
      getData(aderdgCellDescription.getSolution()).data();
  double* aderdgSolutionMin =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
      getData(aderdgCellDescription.getSolutionMin()).data();
  double* aderdgSolutionMax =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
      getData(aderdgCellDescription.getSolutionMax()).data();

  kernels::limiter::generic::c::projectOnADERDGSpace(
      finiteVolumesSolution,
      finiteVolumesSolver->getNumberOfVariables(),
      finiteVolumesSolver->getNodesPerCoordinateAxis(),
      aderdgSolver->getNodesPerCoordinateAxis(),
      aderdgSolution);

  std::fill(aderdgSolutionMin,aderdgSolutionMin+aderdgSolver->getNumberOfVariables(),std::numeric_limits<double>::max());
  std::fill(aderdgSolutionMax,aderdgSolutionMax+aderdgSolver->getNumberOfVariables(),-std::numeric_limits<double>::max()); // !!! Mind the "-"
   // kernels::limiter::generic::c::findCellLocalLimMinAndMax( //TODO Dominic: deprecated method 
      // finiteVolumesSolution,
      // aderdgSolver->getNumberOfVariables(),
      // aderdgSolver->getNodesPerCoordinateAxis(),
      // aderdgSolutionMin,
      // aderdgSolutionMax);
  // This writes the face 0 min max values to arrays for the other faces
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(aderdgSolutionMin, aderdgSolutionMin+aderdgSolver->getNumberOfVariables(),
        aderdgSolutionMin+i*aderdgSolver->getNumberOfVariables());
    std::copy(aderdgSolutionMax, aderdgSolutionMax+aderdgSolver->getNumberOfVariables(),
        aderdgSolutionMax+i*aderdgSolver->getNumberOfVariables());
  }
}

void exahype::solvers::ADERDGAPosterioriSubcellLimiter::couple(
    const int cellDescriptionsIndex,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  exahype::solvers::ADERDGSolver* aderdgSolver =
      static_cast<exahype::solvers::ADERDGSolver*>(exahype::solvers::RegisteredSolvers[_aderdgSolverNumber]);
  exahype::solvers::FiniteVolumesSolver* finiteVolumesSolver =
        static_cast<exahype::solvers::FiniteVolumesSolver*>(exahype::solvers::RegisteredSolvers[_finiteVolumesSolverNumber]);
  int aderdgElement        = aderdgSolver->tryGetElement(cellDescriptionsIndex,_aderdgSolverNumber);
  int finiteVolumesElement = finiteVolumesSolver->tryGetElement(cellDescriptionsIndex,_finiteVolumesSolverNumber);

  if (aderdgElement!=exahype::solvers::Solver::NotFound) {
    auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[aderdgElement];

    assertion2(finiteVolumesElement!=exahype::solvers::Solver::NotFound,cellDescriptionsIndex,_aderdgSolverNumber);
    auto& finiteVolumesCellDescription = exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
        cellDescriptionsIndex)[finiteVolumesElement];

    double* aderdgSolution  = exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolution()).data();
    double* aderdgSolutionMin =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMin()).data();
    double* aderdgSolutionMax =  exahype::solvers::ADERDGSolver::DataHeap::getInstance().
        getData(aderdgCellDescription.getSolutionMax()).data();
    double* finiteVolumesSolution = exahype::solvers::FiniteVolumesSolver::DataHeap::getInstance().
            getData(finiteVolumesCellDescription.getSolution()).data();

    // 1. We first check if the old solution still needs limiting.
    // 1.1 If so, we update the FV solution and project the result back
    // onto the ADER-DG solution space.
    bool oneNeighbourIsTroubledCell = oneNeighbourIsTroubled(aderdgCellDescription.getLimiterStatus());

    if (false) { //TODO Dominic : isTroubled now anticipates the new DG solution
    //if (kernels::limiter::generic::c::isTroubledCell(
    //    aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
    //    aderdgSolutionMin,aderdgSolutionMax)) {
      aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled); // implicit conversion.

      logDebug("couple(...)","Evolving FV solution in cell "<<
          aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<
              " since old solution is still troubled.");

      evolveFiniteVolumesSolution(
          cellDescriptionsIndex,finiteVolumesElement,aderdgElement,
          finiteVolumesSolver,aderdgSolver,
          fineGridVertices,fineGridVerticesEnumerator);
    } else {
      aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok); // implicit conversion.

      if (oneNeighbourIsTroubledCell) {
        logDebug("couple(...)","Evolving FV solution in cell "<<
            aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<
                " since neighbour is troubled.");
        // assertion ( finite volumes solver solution is correct ). todo
        evolveFiniteVolumesSolution(
            cellDescriptionsIndex,finiteVolumesElement,aderdgElement,
            finiteVolumesSolver,aderdgSolver,
            fineGridVertices,fineGridVerticesEnumerator);
      } else {
        // 2. If not so, we first update the solution of the ADER-DG solver, and then
        // check if we need to perform limiting to the new ADER-DG solution.
        // 2.1 If so, we update the FV solution and project the result back
        // onto the ADER-DG solution space.
        logDebug("couple(...)","Evolving ADER-DG solution in cell "<<
            aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<".");

        aderdgSolver->updateSolution(cellDescriptionsIndex,finiteVolumesElement,fineGridVertices,fineGridVerticesEnumerator);

        if(false) { //TODO Dominic: isTroubled now anticipates the new DG solution
        //if (kernels::limiter::generic::c::isTroubledCell(
        //    aderdgSolution,aderdgSolver->getNumberOfVariables(),aderdgSolver->getNodesPerCoordinateAxis(),
        //    aderdgSolutionMin,aderdgSolutionMax)) {
          aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled); // implicit conversion.
          logDebug("couple(...)","Evolved ADER-DG solution in cell " <<
              aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<" needs limiting.");
          evolveFiniteVolumesSolution(
              cellDescriptionsIndex,finiteVolumesElement,aderdgElement,
              finiteVolumesSolver,aderdgSolver,
              fineGridVertices,fineGridVerticesEnumerator);
        } else {
          aderdgCellDescription.setLimiterStatus(exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok); // implicit conversion.

          logDebug("couple(...)","Evolved ADER-DG solution in cell " <<
              aderdgCellDescription.getOffset()<<"x"<<aderdgCellDescription.getSize()<<" does not need limiting.");

          // 3. If everything is fine with the ADER-DG solution, we project it back
          // onto the FV solution space to prepare limiting of the ADER-DG
          // solution values in this cell and in the neighbouring cells.
          kernels::limiter::generic::c::projectOnFVLimiterSpace(
              aderdgSolution,
              aderdgSolver->getNumberOfVariables(),
              aderdgSolver->getNodesPerCoordinateAxis(),
              finiteVolumesSolver->getNodesPerCoordinateAxis(),
              finiteVolumesSolution);
          finiteVolumesCellDescription.setTimeStamp(aderdgCellDescription.getCorrectorTimeStamp());
          finiteVolumesCellDescription.setTimeStepSize(aderdgCellDescription.getCorrectorTimeStepSize());

          // This writes the min max values into an array for face 0;
          std::fill(aderdgSolutionMin,aderdgSolutionMin+aderdgSolver->getNumberOfVariables(),std::numeric_limits<double>::max());
          std::fill(aderdgSolutionMax,aderdgSolutionMax+aderdgSolver->getNumberOfVariables(),-std::numeric_limits<double>::max()); // !!! Mind the "-"
          kernels::limiter::generic::c::findCellLocalMinAndMax(
              aderdgSolution,
              // finiteVolumesSolution, //TODO Dominic: the function doesn't initialize the finiteVolumesSolution anymore
              aderdgSolver->getNumberOfVariables(),
              aderdgSolver->getNodesPerCoordinateAxis(),
              aderdgSolutionMin,
              aderdgSolutionMax);
          // This writes the face 0 min max values to arrays for the other faces
          for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
            std::copy(aderdgSolutionMin, aderdgSolutionMin+aderdgSolver->getNumberOfVariables(),
                aderdgSolutionMin+i*aderdgSolver->getNumberOfVariables());
            std::copy(aderdgSolutionMax, aderdgSolutionMax+aderdgSolver->getNumberOfVariables(),
                aderdgSolutionMax+i*aderdgSolver->getNumberOfVariables());
          }
        }
      }
    }
  }
}
