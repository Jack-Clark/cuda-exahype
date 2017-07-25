/*
 * LimitingADERDGSolver.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#include "LimitingADERDGSolver.h"

#include "tarch/multicore/Lock.h"
#include "kernels/limiter/generic/Limiter.h"

namespace exahype {
namespace solvers {

} /* namespace solvers */
} /* namespace exahype */

tarch::logging::Log exahype::solvers::LimitingADERDGSolver::_log("exahype::solvers::LimitingADERDGSolver");

bool exahype::solvers::LimitingADERDGSolver::limiterDomainOfOneSolverHasChanged() {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG &&
        static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
      return true;
    }
  }
  return false;
}

bool exahype::solvers::LimitingADERDGSolver::isValidCellDescriptionIndex(
    const int cellDescriptionsIndex) const  {
  return _solver->isValidCellDescriptionIndex(cellDescriptionsIndex);
}

exahype::solvers::Solver::SubcellPosition exahype::solvers::LimitingADERDGSolver::computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) {
  return _solver->computeSubcellPositionOfCellOrAncestor(cellDescriptionsIndex,element);
}

exahype::solvers::LimitingADERDGSolver::LimitingADERDGSolver(
    const std::string& identifier,
    std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
    std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
    const double DMPRelaxationParameter,
    const double DMPDifferenceScaling)
    :
    exahype::solvers::Solver(identifier, Solver::Type::LimitingADERDG, solver->getNumberOfVariables(),
        solver->getNumberOfParameters(), solver->getNodesPerCoordinateAxis(), solver->getMaximumMeshSize(),
          solver->getTimeStepping()),
          _solver(std::move(solver)),
          _limiter(std::move(limiter)),
          _limiterDomainHasChanged(false),
          _nextLimiterDomainHasChanged(false),
          _DMPMaximumRelaxationParameter(DMPRelaxationParameter), // TODO externalise
          _DMPDifferenceScaling(DMPDifferenceScaling)           // TODO externalise
{
  assertion(_solver->getNumberOfParameters() == 0);

  assertion(_solver->getTimeStepping()==_limiter->getTimeStepping());
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStamp() const {
  assertionEquals(_limiter->getMinTimeStamp(),_solver->getMinTimeStamp());
  return _solver->getMinTimeStamp();
}

double exahype::solvers::LimitingADERDGSolver::getMinTimeStepSize() const {
  assertionEquals(_limiter->getMinTimeStepSize(),_solver->getMinTimeStepSize());
  return _solver->getMinTimeStepSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinNextTimeStepSize() const {
  assertionEquals(_limiter->getMinNextTimeStepSize(),_solver->getMinNextTimeStepSize());
  return _solver->getMinNextTimeStepSize();
}

void exahype::solvers::LimitingADERDGSolver::updateMinNextTimeStepSize(double value) {
  _solver->updateMinNextTimeStepSize(value);
  _limiter->updateMinNextTimeStepSize(value);
}

void exahype::solvers::LimitingADERDGSolver::initSolverTimeStepData(double value) {
  _solver->initSolverTimeStepData(value);
  _limiter->initSolverTimeStepData(value);
}

void exahype::solvers::LimitingADERDGSolver::synchroniseTimeStepping(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->synchroniseTimeStepping(cellDescriptionsIndex,element);

  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::Troubled
              || solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::NeighbourIsTroubledCell
              || solverPatch.getLimiterStatus()==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);
    _limiter->synchroniseTimeStepping(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::startNewTimeStep() {
  double minNextTimeStepSize =
            std::min( _solver->getMinNextPredictorTimeStepSize(), _limiter->getMinNextTimeStepSize() );

  switch (_timeStepping) {
    case TimeStepping::Global:
      _solver->updateMinNextPredictorTimeStepSize(minNextTimeStepSize);
      _limiter->updateMinNextTimeStepSize(minNextTimeStepSize);

      _solver->startNewTimeStep();
      _limiter->startNewTimeStep();
      break;
    case TimeStepping::GlobalFixed:
      _solver->startNewTimeStep();
      _limiter->startNewTimeStep();
      break;
  } // TODO(Dominic): Switch-case probably not necessary

  _minCellSize     = _nextMinCellSize;
  _maxCellSize     = _nextMaxCellSize;
  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min

  logDebug("startNewTimeStep()","_limiterDomainHasChanged="<<_limiterDomainHasChanged<<",_nextLimiterDomainHasChanged="<<_nextLimiterDomainHasChanged);

  _limiterDomainHasChanged     = _nextLimiterDomainHasChanged;
  _nextLimiterDomainHasChanged = false;
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep() {
  _solver->rollbackToPreviousTimeStep();
  _limiter->rollbackToPreviousTimeStep();

  _nextMinCellSize = std::numeric_limits<double>::max();
  _nextMaxCellSize = -std::numeric_limits<double>::max(); // "-", min
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback() {
  _solver->reconstructStandardTimeSteppingDataAfterRollback();

  _limiter->setMinTimeStamp(_solver->getMinCorrectorTimeStamp());
  _limiter->setMinTimeStepSize(_solver->getMinCorrectorTimeStepSize());
  _limiter->setMinNextTimeStepSize(std::numeric_limits<double>::max());
  _limiter->setPreviousMinTimeStepSize(std::numeric_limits<double>::max());
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseTimeStepData() {
  _solver->reinitialiseTimeStepData();
  _limiter->reinitialiseTimeStepData();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMinCellSize(double minCellSize) {
  _solver->updateNextMinCellSize(minCellSize);
  _limiter->updateNextMinCellSize(minCellSize);

  assertionEquals(_solver->getNextMinCellSize(),_limiter->getNextMinCellSize());
  _nextMinCellSize = _solver->getNextMinCellSize();
}

void exahype::solvers::LimitingADERDGSolver::updateNextMaxCellSize(double maxCellSize) {
  _solver->updateNextMaxCellSize(maxCellSize);
  _limiter->updateNextMaxCellSize(maxCellSize);

  assertionEquals(_solver->getNextMaxCellSize(),_limiter->getNextMaxCellSize());
  _nextMaxCellSize = _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMinCellSize() const {
  assertionEquals(_solver->getNextMinCellSize(),_limiter->getNextMinCellSize());
  return _solver->getNextMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getNextMaxCellSize() const {
  assertionEquals(_solver->getNextMaxCellSize(),_limiter->getNextMaxCellSize());
  return _solver->getNextMaxCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMinCellSize() const {
  assertionEquals(_solver->getMinCellSize(),_limiter->getMinCellSize());
  return _solver->getMinCellSize();
}

double exahype::solvers::LimitingADERDGSolver::getMaxCellSize() const {
  assertionEquals(_solver->getMaxCellSize(),_limiter->getMaxCellSize());
  return _solver->getMaxCellSize();
}

///////////////////////////////////
// MODIFY CELL DESCRIPTION
///////////////////////////////////
bool exahype::solvers::LimitingADERDGSolver::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  return _solver->enterCell(
      fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
      coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
      fineGridPositionOfCell,solverNumber);
}

bool exahype::solvers::LimitingADERDGSolver::leaveCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    const int solverNumber)  {
  return _solver->leaveCell(
          fineGridCell,fineGridVertices,fineGridVerticesEnumerator,
          coarseGridVertices,coarseGridVerticesEnumerator,coarseGridCell,
          fineGridPositionOfCell,solverNumber);
}

///////////////////////////////////
// CELL-LOCAL
//////////////////////////////////
double exahype::solvers::LimitingADERDGSolver::startNewTimeStep(
    const int cellDescriptionsIndex,
    const int element,
    double*   tempEigenvalues)  {
  exahype::solvers::ADERDGSolver::CellDescription& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  double admissibleTimeStepSize = std::numeric_limits<double>::max();

  switch (solverPatch.getLimiterStatus()) {
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok: {
      admissibleTimeStepSize =
          _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      break;
    }
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
      admissibleTimeStepSize =
                  _solver->startNewTimeStep(cellDescriptionsIndex,element,tempEigenvalues);
      int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
      break;
    }
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Troubled:
    case exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell: {
      int limiterElement =
          _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      admissibleTimeStepSize =
          _limiter->startNewTimeStep(cellDescriptionsIndex,limiterElement,tempEigenvalues);
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      LimiterPatch& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      solverPatch.setPreviousCorrectorTimeStepSize(limiterPatch.getPreviousTimeStepSize());
      solverPatch.setCorrectorTimeStamp(limiterPatch.getTimeStamp());
      solverPatch.setCorrectorTimeStepSize(limiterPatch.getTimeStepSize());

      solverPatch.setPredictorTimeStamp(limiterPatch.getTimeStamp()+limiterPatch.getTimeStepSize());
      solverPatch.setPredictorTimeStepSize(limiterPatch.getTimeStepSize()); // TODO(Dominic): Reassess this for Standard time stepping.
      break;
    }
  }

  return admissibleTimeStepSize;
}

void exahype::solvers::LimitingADERDGSolver::rollbackToPreviousTimeStep(
    const int cellDescriptionsIndex,
    const int solverElement) {
  _solver->rollbackToPreviousTimeStep(cellDescriptionsIndex,solverElement);

  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLimiterStatus()==SolverPatch::Troubled
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsTroubledCell
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsNeighbourOfTroubledCell) {
    int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    _limiter->rollbackToPreviousTimeStep(cellDescriptionsIndex,limiterElement);
  }
}

void exahype::solvers::LimitingADERDGSolver::reconstructStandardTimeSteppingDataAfterRollback(
    const int cellDescriptionsIndex,
    const int solverElement) const {
  _solver->reconstructStandardTimeSteppingDataAfterRollback(cellDescriptionsIndex,solverElement);

  SolverPatch& solverPatch =
        _solver->getCellDescription(cellDescriptionsIndex,solverElement);
  if (solverPatch.getLimiterStatus()==SolverPatch::Troubled
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsTroubledCell
      || solverPatch.getLimiterStatus()==SolverPatch::NeighbourIsNeighbourOfTroubledCell) {
    int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    LimiterPatch& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    limiterPatch.setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
    limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
    limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
  }
}

void exahype::solvers::LimitingADERDGSolver::setInitialConditions(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  _solver->setInitialConditions(
      cellDescriptionsIndex,element,
      fineGridVertices,fineGridVerticesEnumerator);
}

void exahype::solvers::LimitingADERDGSolver::initialiseLimiter( // TODO(Dominic): Remove is not used anymore
    const int cellDescriptionsIndex,
    const int element,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  int limiterElement         = exahype::solvers::Solver::NotFound;
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    // do nothing
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    // TODO(Dominic): Use solver patch's cell type
    // TODO(Dominic): This is some sort of mesh refinement.
    // In AMR settings, we need to add ancestor and descendant cells
    // after we add a limiter cell description (this will be fun).
    fineGridCell.addNewCellDescription(
        solverPatch.getSolverNumber(),
        LimiterPatch::Type::Cell,
        solverPatch.getLevel(),
        solverPatch.getParentIndex(),
        solverPatch.getSize(),
        solverPatch.getOffset());

    limiterElement = tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    limiterPatch =
        &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
    _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);

    limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
    _limiter->setInitialConditions(cellDescriptionsIndex,limiterElement,fineGridVertices,fineGridVerticesEnumerator);

    // Finite Volumes -> ADER-DG
    solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
    kernels::limiter::generic::c::projectOnDGSpace(
       limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    // TODO(Dominic): Use solver patch's cell type
    // TODO(Dominic): This is some sort of mesh refinement.
    // In AMR settings, we need to add ancestor and descendant cells
    // after we add a limiter cell description (this will be fun).
    fineGridCell.addNewCellDescription(
        solverPatch.getSolverNumber(),
        LimiterPatch::Type::Cell,
        solverPatch.getLevel(),
        solverPatch.getParentIndex(),
        solverPatch.getSize(),
        solverPatch.getOffset());

    limiterElement = tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    limiterPatch =
        &LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex())[limiterElement];
    _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);
    limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();

    // ADER-DG -> Finite Volumes
    solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        limiterSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    break;
  }
}

/**
 * This method assumes the ADERDG solver's cell-local limiter status has
 * already been determined.
 */
void exahype::solvers::LimitingADERDGSolver::updateSolution(
    const int cellDescriptionsIndex,
    const int element,
    double** tempStateSizedVectors,
    double** tempUnknowns,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator)  {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  // 1. Update the solution in the cells
  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    break;
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
    _solver->updateSolution(
        cellDescriptionsIndex,element,
        tempStateSizedVectors,tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();

    _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);
    limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnFVLimiterSpace(
        solverSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        limiterSolution);
    } break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
    _limiter->updateSolution(
        cellDescriptionsIndex,limiterElement,
        tempStateSizedVectors,
        tempUnknowns,
        fineGridVertices,fineGridVerticesEnumerator);
    assertion(limiterElement!=exahype::solvers::Solver::NotFound);
    exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
        _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

    solverSolution = DataHeap::getInstance().getData(
        solverPatch.getSolution()).data();
    limiterSolution = DataHeap::getInstance().getData(
        limiterPatch.getSolution()).data();

    // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
    kernels::limiter::generic::c::projectOnDGSpace(limiterSolution,_solver->getNumberOfVariables(),
        _solver->getNodesPerCoordinateAxis(),
        _limiter->getGhostLayerWidth(),
        solverSolution);
    } break;
  }
}

//void printSolutionMinOrMax(const double* minOrMax,const int numberOfVariables,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int i = 0; i < DIMENSIONS_TIMES_TWO*numberOfVariables; ++i) {
//    std::cout << minOrMax[i] << ",";
//  }
//  std::cout << std::endl;
//}
//
//void printNormalFluxes(const double* flux,const int numberOfFaceUnknowns,const char* identifier) {
//  std::cout << identifier << "=" ;
//
//  for (int d=0; d<DIMENSIONS_TIMES_TWO;d++) {
//    for (int i = 0; i < numberOfFaceUnknowns; ++i) {
//      std::cout << flux[i+d*numberOfFaceUnknowns] << ",";
//    }
//    std::cout << "|";
//  }
//  std::cout << std::endl;
//}

bool exahype::solvers::LimitingADERDGSolver::updateMergedLimiterStatusAndMinAndMaxAfterSolutionUpdate(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  bool solutionIsValid = evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(solverPatch) &&
                         evaluatePhysicalAdmissibilityCriterion(solverPatch); // after min and max was found


  switch (solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
        const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
        assertion2(limiterElement!=exahype::solvers::Solver::NotFound,limiterElement,cellDescriptionsIndex);
        LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
        determineLimiterMinAndMax(solverPatch,limiterPatch);
      }
      break;
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      // Already computed.
      break;
  }

  return determineMergedLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid);
}

bool exahype::solvers::LimitingADERDGSolver::updateMergedLimiterStatusAndMinAndMaxAfterSetInitialConditions(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  determineSolverMinAndMax(solverPatch);

  bool solutionIsValid = evaluatePhysicalAdmissibilityCriterion(solverPatch); // Only evaluate the PAD

  return determineMergedLimiterStatusAfterSolutionUpdate(solverPatch,!solutionIsValid);
}

bool exahype::solvers::LimitingADERDGSolver::determineMergedLimiterStatusAfterSolutionUpdate(
    SolverPatch& solverPatch,const bool isTroubled) const {
  bool limiterDomainHasChanged=false;

  switch (solverPatch.getLimiterStatus()) {
  case SolverPatch::LimiterStatus::Ok:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    if (isTroubled) {
      limiterDomainHasChanged = true;

      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
      }
    } else {
      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
      }
    }
    break;
  case SolverPatch::LimiterStatus::Troubled:
    if (!isTroubled) {
      limiterDomainHasChanged = true;

      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
      }
    } else {
      for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
        solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
      }
    }
    break;
  default:
    break;
  }

//  // TODO(Dominic): Remove
//  tarch::la::Vector<DIMENSIONS,double> poi1(0.462963,0.833333);
//  tarch::la::Vector<DIMENSIONS,double> poi2(0.364198,0.808642);
//  if (tarch::la::equals(solverPatch.getOffset(),poi1,1e-3)
//  || tarch::la::equals(solverPatch.getOffset(),poi2,1e-3)
//  ) {
//    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
//      solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Troubled);
//    }
//  } else {
//    for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
//      solverPatch.setMergedLimiterStatus(i,SolverPatch::LimiterStatus::Ok);
//    }
//  }
//
//  // TODO(Dominic): Remove
//  limiterDomainHasChanged = true;

  return limiterDomainHasChanged;
}


bool exahype::solvers::LimitingADERDGSolver::evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // 1. Check if the DMP is satisfied and search for the min and max
  // Write the new min and max to the storage reserved for face 0
  bool dmpIsSatisfied = kernels::limiter::generic::c::discreteMaximumPrincipleAndMinAndMaxSearch(
        solution,_solver->getNumberOfVariables(),_solver->getNodesPerCoordinateAxis(),
        _DMPMaximumRelaxationParameter, _DMPDifferenceScaling,
        solutionMin,solutionMax);

  // 2. Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(
        solutionMin,
        solutionMin+_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy(
        solutionMax,
        solutionMax+_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
  }

  if (dmpIsSatisfied) {
    return true;
  }

  return false;
}

bool exahype::solvers::LimitingADERDGSolver::evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch) {
  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  return _solver->physicalAdmissibilityDetection(solutionMin,solutionMax);
}

void exahype::solvers::LimitingADERDGSolver::determineMinAndMax(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch =
      _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      determineSolverMinAndMax(solverPatch);
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion1(limiterElement!=exahype::solvers::Solver::NotFound,solverPatch.toString());
      exahype::solvers::FiniteVolumesSolver::CellDescription& limiterPatch =
          _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      determineLimiterMinAndMax(solverPatch,limiterPatch);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::determineSolverMinAndMax(SolverPatch& solverPatch) {
  double* solution = DataHeap::getInstance().getData(
      solverPatch.getSolution()).data();

  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalMinAndMax(
      solution,_numberOfVariables,_nodesPerCoordinateAxis,solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(
        solutionMin,
        solutionMin+_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy(
        solutionMax,
        solutionMax+_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
  }
}

void exahype::solvers::LimitingADERDGSolver::determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch) {
  double* limiterSolution = DataHeap::getInstance().getData(
      limiterPatch.getSolution()).data();

  double* solutionMin = DataHeap::getInstance().getData(
      solverPatch.getSolutionMin()).data();
  double* solutionMax = DataHeap::getInstance().getData(
      solverPatch.getSolutionMax()).data();

  // Write the result to the face with index "0"
  kernels::limiter::generic::c::findCellLocalLimiterMinAndMax(
      limiterSolution,_numberOfVariables,_nodesPerCoordinateAxis,_limiter->getGhostLayerWidth(),solutionMin,solutionMax);

  // Copy the result on the other faces as well
  for (int i=1; i<DIMENSIONS_TIMES_TWO; ++i) {
    std::copy(
        solutionMin,
        solutionMin+_numberOfVariables, // past-the-end element
        solutionMin+i*_numberOfVariables);
    std::copy(
        solutionMax,
        solutionMax+_numberOfVariables, // past-the-end element
        solutionMax+i*_numberOfVariables);
  }
}

/**
 * Determine the limiter status of a cell description after a solution update.
 * Check which cell is troubled and which cell holds a valid solution.
 * If the limiter subdomain changes, i.e., if a cell changes from holding a
 * valid solution (Ok) to troubled (Troubled) or vice versa,
 * this function returns true.
 *
 * Further reinitialise the merged limiter status on every
 * face with either Ok or Troubled.
 *
 * Do not override the cell-based limiter status of the solver only
 * the merged ones on the faces. We need to hold this value a
 * little longer.
 */
void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    SolverPatch& solverPatch,
    const int faceIndex,
    const SolverPatch::LimiterStatus& neighbourLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      switch (neighbourLimiterStatus) {
        case SolverPatch::LimiterStatus::Troubled:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
          break;
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (neighbourLimiterStatus) {
        case SolverPatch::LimiterStatus::Troubled:
          solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourLimiterStatus(
    SolverPatch& solverPatch,
    const int faceIndex,
    const SolverPatch::LimiterStatus& neighbourLimiterStatus,
    const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const {
  SolverPatch::LimiterStatus limiterStatus = determineLimiterStatus(solverPatch);

  mergeWithNeighbourLimiterStatus(solverPatch,faceIndex,neighbourLimiterStatus);

  if (limiterStatus==SolverPatch::LimiterStatus::Ok) {
    if (neighbourLimiterStatus==SolverPatch::LimiterStatus::Ok
        && neighbourOfNeighbourLimiterStatus==SolverPatch::LimiterStatus::Troubled) {
      solverPatch.setMergedLimiterStatus(faceIndex,SolverPatch::LimiterStatus::NeighbourIsTroubledCell);
    }
  }
}

/**
 * Iterate over the merged limiter statuses per face and
 * determine a unique value.
 */
exahype::solvers::LimitingADERDGSolver::SolverPatch::LimiterStatus
exahype::solvers::LimitingADERDGSolver::determineLimiterStatus(
    SolverPatch& solverPatch) const {
  // 1. Determine new limiter status.
  SolverPatch::LimiterStatus limiterStatus = SolverPatch::LimiterStatus::Ok;
  for (int i=0; i<DIMENSIONS_TIMES_TWO; i++) {
    switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      switch (solverPatch.getMergedLimiterStatus(i)) {
      case SolverPatch::LimiterStatus::Troubled:
        limiterStatus = SolverPatch::LimiterStatus::Troubled;
        break;
      case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsTroubledCell;
        break;
      case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell;
        break;
      default:
        break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (solverPatch.getMergedLimiterStatus(i)) {
      case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        limiterStatus = SolverPatch::LimiterStatus::NeighbourIsTroubledCell;
        break;
      default:
        break;
      }
      break;
      default:
        break;
    }
  }
  return limiterStatus;
}

void exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers(
    const int cellDescriptionsIndex,
    const int element,
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus         = determineLimiterStatus(solverPatch);

  int limiterElement = _limiter->tryGetElement(
          fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
  case SolverPatch::LimiterStatus::Ok:
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      assertion(previousLimiterStatus==SolverPatch::LimiterStatus::Troubled ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsTroubledCell ||
          previousLimiterStatus==SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell);

      limiterPatch = &_limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);
      limiterPatch->setType(LimiterPatch::Type::Erased);
      _limiter->ensureNoUnnecessaryMemoryIsAllocated(*limiterPatch);
      LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).erase(
          LimiterHeap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).begin()+limiterElement);
    }
    break;
  case SolverPatch::LimiterStatus::Troubled:
  case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
  case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
    // TODO(Dominic): Add to docu: We need to rollback the solution for this case since we need to supply the
    // NeighbourIsTroubledCell neighbours with solution values from the old time step.
    switch (previousLimiterStatus) {
    case SolverPatch::LimiterStatus::Troubled:  // TODO(Dominic): Add to docu: Here we work with a valid old FVM solution
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->rollbackSolution(
          fineGridCell.getCellDescriptionsIndex(),limiterElement,
          fineGridVertices,fineGridVerticesEnumerator);
      break;
    case SolverPatch::LimiterStatus::Ok: // TODO(Dominic): Add to docu: Here we work with a valid old ADER-DG solution
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: // TODO(Dominic): Possibly dangerous if nans appear
      _solver->rollbackSolution(
                fineGridCell.getCellDescriptionsIndex(),element,
                fineGridVertices,fineGridVerticesEnumerator);

      if (limiterElement==exahype::solvers::Solver::NotFound) {
        assertion(previousLimiterStatus==SolverPatch::LimiterStatus::Ok);

        // TODO(Dominic): Use solver patch's cell type
        // TODO(Dominic): This is some sort of mesh refinement.
        // In AMR settings, we need to add ancestor and descendant cells
        // after we add a limiter cell description (this will be fun).
        fineGridCell.addNewCellDescription(
            solverPatch.getSolverNumber(),
            LimiterPatch::Type::Cell,
            solverPatch.getLevel(),
            solverPatch.getParentIndex(),
            solverPatch.getSize(),
            solverPatch.getOffset());

        limiterElement = _limiter->tryGetElement(
            fineGridCell.getCellDescriptionsIndex(),solverPatch.getSolverNumber());

        solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

        limiterPatch = &_limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);
        _limiter->ensureNecessaryMemoryIsAllocated(*limiterPatch);

        limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
        kernels::limiter::generic::c::projectOnFVLimiterSpace(
            solverSolution,_solver->getNumberOfVariables(),
            _solver->getNodesPerCoordinateAxis(),
            _limiter->getGhostLayerWidth(),
            limiterSolution); // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      } else {
        limiterPatch = &_limiter->getCellDescription(fineGridCell.getCellDescriptionsIndex(),limiterElement);
        _limiter->rollbackSolution(cellDescriptionsIndex,limiterElement,fineGridVertices,fineGridVerticesEnumerator);
      }

      // Copy time step data from the solver patch
      limiterPatch->setPreviousTimeStepSize(solverPatch.getPreviousCorrectorTimeStepSize());
      limiterPatch->setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch->setTimeStepSize(solverPatch.getCorrectorTimeStepSize());

      // Copy more geometry information from the solver patch
      limiterPatch->setIsInside(solverPatch.getIsInside());
      break;
    }
    break;
  }
}

void exahype::solvers::LimitingADERDGSolver::recomputeSolution(
        const int cellDescriptionsIndex, const int element,
         double** tempStateSizedArrays,
         double** tempUnknowns,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus previousLimiterStatus = solverPatch.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus         = determineLimiterStatus(solverPatch);

  int limiterElement =
      _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  LimiterPatch* limiterPatch = nullptr;

  double* solverSolution  = nullptr;
  double* limiterSolution = nullptr;

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
      assertion(limiterElement==exahype::solvers::Solver::NotFound);
      // do nothing
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: // TODO(Dominic): Add to docu: Here, were just went back one time step to supply the NT neighbours with old limiter unknowns.
      switch (previousLimiterStatus) {
        case SolverPatch::LimiterStatus::Ok:
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          // TODO(Dominic): Old code. Keep a while for reference.
          //          _solver->addUpdateToSolution(solverPatch,fineGridVertices,fineGridVerticesEnumerator); // It must always be addUpdateToSolution here since we do not want to use compute the surface integral again
          _solver->swapSolutionAndPreviousSolution(cellDescriptionsIndex,element);

          assertion(limiterElement!=exahype::solvers::Solver::NotFound);
          limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
          solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();
          limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();

          // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
          kernels::limiter::generic::c::projectOnFVLimiterSpace(
              solverSolution,_solver->getNumberOfVariables(),
              _solver->getNodesPerCoordinateAxis(),
              _limiter->getGhostLayerWidth(),
              limiterSolution);
          break;
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
        case SolverPatch::LimiterStatus::Troubled:
          _limiter->swapSolutionAndPreviousSolution(cellDescriptionsIndex,limiterElement);

          limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
          limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
          solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

          // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
          kernels::limiter::generic::c::projectOnDGSpace(
              limiterSolution,_solver->getNumberOfVariables(),
              _solver->getNodesPerCoordinateAxis(),
              _limiter->getGhostLayerWidth(),
              solverSolution);
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);

      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      assertionEquals1(limiterPatch->getPreviousTimeStepSize(),solverPatch.getPreviousCorrectorTimeStepSize(),cellDescriptionsIndex);
      assertionEquals1(limiterPatch->getTimeStamp(),solverPatch.getCorrectorTimeStamp(),cellDescriptionsIndex);
      assertionEquals1(limiterPatch->getTimeStepSize(),solverPatch.getCorrectorTimeStepSize(),cellDescriptionsIndex);

      // 1. Evolve solution to desired  time step again
      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      // 2. Project FV solution on ADER-DG space
      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);

      limiterPatch = &_limiter->getCellDescription(cellDescriptionsIndex,limiterElement);

      _limiter->updateSolution(cellDescriptionsIndex,limiterElement,
                               tempStateSizedArrays,tempUnknowns,
                               fineGridVertices,fineGridVerticesEnumerator);

      limiterSolution = DataHeap::getInstance().getData(limiterPatch->getSolution()).data();
      solverSolution  = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

      // TODO(Dominic): Add virtual method. The current implementation depends on a particular kernel.
      kernels::limiter::generic::c::projectOnDGSpace(
          limiterSolution,_solver->getNumberOfVariables(),
          _solver->getNodesPerCoordinateAxis(),
          _limiter->getGhostLayerWidth(),
          solverSolution);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::prepareNextNeighbourMerging(
    const int cellDescriptionsIndex,const int element,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const {
  _solver->prepareNextNeighbourMerging(
      cellDescriptionsIndex,element,
      fineGridVertices,fineGridVerticesEnumerator);

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->prepareNextNeighbourMerging(
        cellDescriptionsIndex,limiterElement,
        fineGridVertices,fineGridVerticesEnumerator);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateMergedLimiterStatus(
    const int cellDescriptionsIndex,const int element) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus mergedLimiterStatus = determineLimiterStatus(solverPatch);

  for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
    solverPatch.setMergedLimiterStatus(i,mergedLimiterStatus);
  }
}

void exahype::solvers::LimitingADERDGSolver::updateLimiterStatus(
    const int cellDescriptionsIndex,const int element) const  {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  solverPatch.setLimiterStatus(determineLimiterStatus(solverPatch));
}

void exahype::solvers::LimitingADERDGSolver::preProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->preProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->preProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::postProcess(
    const int cellDescriptionsIndex,
    const int element) {
  _solver->postProcess(cellDescriptionsIndex,element);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  _limiter->postProcess(cellDescriptionsIndex,limiterElement);
}

void exahype::solvers::LimitingADERDGSolver::prolongateDataAndPrepareDataRestriction(
    const int cellDescriptionsIndex,
    const int element) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = exahype::solvers::Solver::NotFound;

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,element);

      limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,limiterElement);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->prolongateDataAndPrepareDataRestriction(cellDescriptionsIndex,limiterElement);
      break;
  }
}

void exahype::solvers::LimitingADERDGSolver::restrictData(
    const int cellDescriptionsIndex,
    const int element,
    const int parentCellDescriptionsIndex,
    const int parentElement,
    const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement       = exahype::solvers::Solver::NotFound;
  int parentLimiterElement = exahype::solvers::Solver::NotFound;

  switch(solverPatch.getLimiterStatus()) {
    case SolverPatch::LimiterStatus::Ok:
      _solver->restrictData(cellDescriptionsIndex,element,parentCellDescriptionsIndex,parentElement,subcellIndex);
      break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      _solver->restrictData(cellDescriptionsIndex,element,parentCellDescriptionsIndex,parentElement,subcellIndex);

      limiterElement       = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      parentLimiterElement = tryGetLimiterElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->restrictData(cellDescriptionsIndex,limiterElement,parentCellDescriptionsIndex,parentLimiterElement,subcellIndex);
      break;
    case SolverPatch::LimiterStatus::Troubled:
      limiterElement       = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      parentLimiterElement = tryGetLimiterElement(parentCellDescriptionsIndex,solverPatch.getSolverNumber());
      _limiter->restrictData(cellDescriptionsIndex,limiterElement,parentCellDescriptionsIndex,parentLimiterElement,subcellIndex);
      break;
  }
}


///////////////////////////////////
// NEIGHBOUR
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeLimiterStatusOfNeighbours(
      const int                                 SolverPatchsIndex1,
      const int                                 element1,
      const int                                 SolverPatchsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int SolverPatchsIndexLeft  = SolverPatchsIndex1;
  int elementLeft            = element1;
  int faceIndexLeft          = faceIndex1;

  int SolverPatchsIndexRight = SolverPatchsIndex2;
  int elementRight           = element2;
  int faceIndexRight         = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    SolverPatchsIndexLeft  = SolverPatchsIndex2;
    elementLeft            = element2;
    faceIndexLeft          = faceIndex2;

    SolverPatchsIndexRight = SolverPatchsIndex1;
    elementRight           = element1;
    faceIndexRight         = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(SolverPatchsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(SolverPatchsIndexRight,elementRight);

  // We need to copy the limiter status since the routines below modify
  // the limiter status on the cell descriptions.
  const SolverPatch::LimiterStatus& limiterStatusLeft  = solverPatchLeft.getMergedLimiterStatus(faceIndexLeft);
  const SolverPatch::LimiterStatus& limiterStatusRight = solverPatchRight.getMergedLimiterStatus(faceIndexRight);
  mergeWithNeighbourLimiterStatus(solverPatchLeft,faceIndexLeft,limiterStatusRight);
  mergeWithNeighbourLimiterStatus(solverPatchRight,faceIndexRight,limiterStatusLeft);
}

// TODO(Dominic): Remove
//void exahype::solvers::LimitingADERDGSolver::mergeLimiterStatusOfNeighboursOfNeighbours(
//      const int                                 SolverPatchsIndex1,
//      const int                                 element1,
//      const int                                 SolverPatchsIndex2,
//      const int                                 element2,
//      const tarch::la::Vector<DIMENSIONS, int>& pos1,
//      const tarch::la::Vector<DIMENSIONS, int>& pos2) {
//  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));
//
//  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
//  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
//  const int faceIndex1 = 2 * normalDirection +
//      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
//  const int faceIndex2 = 2 * normalDirection +
//      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!
//
//  int SolverPatchsIndexLeft  = SolverPatchsIndex1;
//  int elementLeft            = element1;
//  int faceIndexLeft          = faceIndex1;
//
//  int SolverPatchsIndexRight = SolverPatchsIndex2;
//  int elementRight           = element2;
//  int faceIndexRight         = faceIndex2;
//
//  if (pos1(normalDirection) > pos2(normalDirection)) {
//    SolverPatchsIndexLeft  = SolverPatchsIndex2;
//    elementLeft            = element2;
//    faceIndexLeft          = faceIndex2;
//
//    SolverPatchsIndexRight = SolverPatchsIndex1;
//    elementRight           = element1;
//    faceIndexRight         = faceIndex1;
//  }
//
//  SolverPatch& solverPatchLeft  = _solver->getCellDescription(SolverPatchsIndexLeft,elementLeft);
//  SolverPatch& solverPatchRight = _solver->getCellDescription(SolverPatchsIndexRight,elementRight);
//
//  // We need to copy the limiter status since the routines below modify
//  // the limiter status on the cell descriptions.
//  const SolverPatch::LimiterStatus& limiterStatusLeft       = solverPatchLeft.getMergedLimiterStatus(faceIndexLeft);
//  const SolverPatch::LimiterStatus& limiterStatusLeftLeft   = solverPatchLeft.getMergedLimiterStatus(faceIndexRight);
//  const SolverPatch::LimiterStatus& limiterStatusRight      = solverPatchRight.getMergedLimiterStatus(faceIndexRight);
//  const SolverPatch::LimiterStatus& limiterStatusRightRight = solverPatchRight.getMergedLimiterStatus(faceIndexLeft);
//  mergeWithNeighbourLimiterStatus(solverPatchLeft,faceIndexLeft,limiterStatusRight, limiterStatusRightRight);
//  mergeWithNeighbourLimiterStatus(solverPatchRight,faceIndexRight,limiterStatusLeft, limiterStatusLeftLeft);
//}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch& pLeft,
  SolverPatch& pRight,
  const int faceIndexLeft,
  const int faceIndexRight
) const {
  if (pLeft.getType()==SolverPatch::Cell ||
      pRight.getType()==SolverPatch::Cell) {
    assertion( pLeft.getSolverNumber() == pRight.getSolverNumber() );
    const int numberOfVariables = getNumberOfVariables();
    double* minLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMin()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* minRight = DataHeap::getInstance().getData( pRight.getSolutionMin()  ).data() + faceIndexRight * numberOfVariables;
    double* maxLeft  = DataHeap::getInstance().getData( pLeft.getSolutionMax()  ).data()  + faceIndexLeft  * numberOfVariables;
    double* maxRight = DataHeap::getInstance().getData( pRight.getSolutionMax()  ).data() + faceIndexRight * numberOfVariables;

    for (int i=0; i<numberOfVariables; i++) {
      const double min = std::min(
          *(minLeft+i),
          *(minRight+i)
      );
      const double max = std::max(
          *(maxLeft+i),
          *(maxRight+i)
      );

      *(minLeft+i)  = min;
      *(minRight+i) = min;

      *(maxLeft+i)  = max;
      *(maxRight+i) = max;
    }
  } // else do nothing
}

void exahype::solvers::LimitingADERDGSolver::mergeNeighbours(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1. Communicate the data between the solvers that is necessary for the solves
  mergeNeighboursBasedOnLimiterStatus(
      cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,
      pos1,pos2,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);

  // 2. Merge the min and max of both cell description's solver's
  // solution value.
  // 2.1 Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 2.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

// TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
// They depend on the isRecomputation value
void exahype::solvers::LimitingADERDGSolver::mergeNeighboursBasedOnLimiterStatus(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    const bool                                isRecomputation,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  // 1. Merge solver solution or limiter solution values in
  // non-overlapping parts of solver and limiter domain:
  SolverPatch& solverPatch1 = _solver->getCellDescription(cellDescriptionsIndex1,element1);
  SolverPatch& solverPatch2 = _solver->getCellDescription(cellDescriptionsIndex2,element2);
  int limiterElement1 = tryGetLimiterElement(cellDescriptionsIndex1,solverPatch1.getSolverNumber());
  int limiterElement2 = tryGetLimiterElement(cellDescriptionsIndex2,solverPatch2.getSolverNumber());

  SolverPatch::LimiterStatus limiterStatus1 = solverPatch1.getLimiterStatus();
  SolverPatch::LimiterStatus limiterStatus2 = solverPatch2.getLimiterStatus();

  if (isRecomputation) {
    limiterStatus1 = solverPatch1.getMergedLimiterStatus(0);
    limiterStatus2 = solverPatch2.getMergedLimiterStatus(0);

    for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
      assertionEquals(limiterStatus1,solverPatch1.getMergedLimiterStatus(i));
      assertionEquals(limiterStatus2,solverPatch2.getMergedLimiterStatus(i));
    } // Dead code elimination will get rid of this loop if Asserts flag is not set.
  }

  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Ok:
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::Troubled:
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          break;
        default:
          break;
      }
      break;
  }
  // 2. Merge limiter solution values in overlapping part
  // of solver and limiter domain:
  switch (limiterStatus1) {
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      switch (limiterStatus2) {
        case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
          _limiter->mergeNeighbours(cellDescriptionsIndex1,limiterElement1,cellDescriptionsIndex2,limiterElement2,pos1,pos2,
                                    tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          if (!isRecomputation) {
            _solver->mergeNeighbours(cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices); // Which one is left and right is checked internally again.
          }
          break;
        default:
          break;
      }
      break;
    default:
      break;
  }
  // !!! Note that it is not necessary to update the boundary layer of the solver patches with status NeighbourIsNeighbourOfTroubledCell however
  // we want to use the symmetric mergeNeighbours method of the finite volumes solver.
  // The non-overlapping part between limiter and solver domains is thus skipped in the parallel mergeWithNeighbourData methods.
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMax(
    const int                                 cellDescriptionsIndex1,
    const int                                 element1,
    const int                                 cellDescriptionsIndex2,
    const int                                 element2,
    const tarch::la::Vector<DIMENSIONS, int>& pos1,
    const tarch::la::Vector<DIMENSIONS, int>& pos2,
    double**                                  tempFaceUnknowns,
    double**                                  tempStateSizedVectors,
    double**                                  tempStateSizedSquareMatrices) {
  assertion1(tarch::la::countEqualEntries(pos1,pos2)==(DIMENSIONS-1),tarch::la::countEqualEntries(pos1,pos2));

  // 1.1. Determine "left" and "right" solverPatch
  const int normalDirection = tarch::la::equalsReturnIndex(pos1, pos2);
  assertion(normalDirection >= 0 && normalDirection < DIMENSIONS);
  const int faceIndex1 = 2 * normalDirection +
      (pos2(normalDirection) > pos1(normalDirection) ? 1 : 0); // !!! Be aware of the ">" !!!
  const int faceIndex2 = 2 * normalDirection +
      (pos1(normalDirection) > pos2(normalDirection) ? 1 : 0);   // !!! Be aware of the ">" !!!

  int cellDescriptionsIndexLeft  = cellDescriptionsIndex1;
  int elementLeft                = element1;
  int faceIndexLeft              = faceIndex1;

  int cellDescriptionsIndexRight = cellDescriptionsIndex2;
  int elementRight               = element2;
  int faceIndexRight             = faceIndex2;

  if (pos1(normalDirection) > pos2(normalDirection)) {
    cellDescriptionsIndexLeft  = cellDescriptionsIndex2;
    elementLeft                = element2;
    faceIndexLeft              = faceIndex2;

    cellDescriptionsIndexRight = cellDescriptionsIndex1;
    elementRight               = element1;
    faceIndexRight             = faceIndex1;
  }

  SolverPatch& solverPatchLeft  = _solver->getCellDescription(cellDescriptionsIndexLeft,elementLeft);
  SolverPatch& solverPatchRight = _solver->getCellDescription(cellDescriptionsIndexRight,elementRight);

  // 1.2. Merge min/max of both solver patches
  mergeSolutionMinMaxOnFace(solverPatchLeft,solverPatchRight,faceIndexLeft,faceIndexRight);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);

  mergeWithBoundaryDataBasedOnLimiterStatus(
      cellDescriptionsIndex,element,solverPatch.getLimiterStatus(),posCell,posBoundary,
      false, /* isRecomputation */
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithBoundaryDataBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const SolverPatch::LimiterStatus&         limiterStatus,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell:
      assertion(limiterStatus==SolverPatch::LimiterStatus::Ok || limiterElement!=exahype::solvers::Solver::NotFound);
      if (!isRecomputation) {
        _solver->mergeWithBoundaryData(cellDescriptionsIndex,element,posCell,posBoundary,
                                       tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      }
      break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell:
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithBoundaryData(cellDescriptionsIndex,limiterElement,posCell,posBoundary,
                                     tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices);
      break;
    default:
      break;
  }
}

#ifdef Parallel
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerNeighbourCommunication    = 1;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerForkOrJoinCommunication   = 0;
const int exahype::solvers::LimitingADERDGSolver::DataMessagesPerMasterWorkerCommunication = 0;

///////////////////////////////////
// NEIGHBOUR - Time marching
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMetadata(
      const int neighbourTypeAsInt,
      const int cellDescriptionsIndex,
      const int element) {
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithNeighbourMetadata(neighbourTypeAsInt,cellDescriptionsIndex,limiterElement);
  } // There is no drop method for metadata necessary

  _solver->mergeWithNeighbourMetadata(neighbourTypeAsInt,cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  sendMinAndMaxToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);

  sendDataToNeighbourBasedOnLimiterStatus(
      toRank,cellDescriptionsIndex,element,src,dest,
      false,/* isRecomputation */
      x,level);
}

void exahype::solvers::LimitingADERDGSolver::sendMinAndMaxToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
                                                                          // |src|dest| : 1; |dest|src| : 0

  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMin()));
  assertion(DataHeap::getInstance().isValidIndex(solverPatch.getSolutionMax()));

  // We append all the max values to the min values.
  // And then append the limiter status as double
  std::vector<double> minAndMaxToSend(2*_numberOfVariables);
  for (int i=0; i<_numberOfVariables; i++) {
    minAndMaxToSend[i]                    = DataHeap::getInstance().getData( solverPatch.getSolutionMin() )[faceIndex*_numberOfVariables+i];
    minAndMaxToSend[i+_numberOfVariables] = DataHeap::getInstance().getData( solverPatch.getSolutionMax() )[faceIndex*_numberOfVariables+i];
  }

  DataHeap::getInstance().sendData(
      minAndMaxToSend, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}


void exahype::solvers::LimitingADERDGSolver::sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();
  if (isRecomputation) {
    limiterStatus = solverPatch.getMergedLimiterStatus(0);

    for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
      assertionEquals(limiterStatus,solverPatch.getMergedLimiterStatus(i));
    } // Dead code elimination will get rid of this loop if Asserts flag is not set.
  }

  logDebug("sendDataToNeighbourBasedOnLimiterStatus(...)", "send data for solver " << _identifier << " from rank " <<
               toRank << " at vertex x=" << x << ", level=" << level <<
               ", src=" << src << ", dest=" << dest <<", limiterStatus="<<limiterStatus);

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok: {
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      } break;
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::Troubled: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _solver->sendDataToNeighbour(toRank,cellDescriptionsIndex,element,src,dest,x,level);
      _limiter->sendDataToNeighbour(toRank,cellDescriptionsIndex,limiterElement,src,dest,x,level);
      } break;
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  std::vector<double> emptyMessage(0);
  for(int sends=0; sends<DataMessagesPerNeighbourCommunication; ++sends)
    DataHeap::getInstance().sendData(
        emptyMessage, toRank, x, level,
        peano::heap::MessageType::NeighbourCommunication);

  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourData(
    const int                                    fromRank,
    const int                                    neighbourTypeAsInt,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  mergeWithNeighbourDataBasedOnLimiterStatus(
      fromRank,neighbourTypeAsInt,cellDescriptionsIndex,element,src,dest,
      false,/*isRecomputation*/
      tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);

  mergeWithNeighbourMinAndMax(fromRank,cellDescriptionsIndex,element,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourDataBasedOnLimiterStatus(
    const int                                    fromRank,
    const int                                    neighbourTypeAsInt,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const bool                                   isRecomputation,
    double**                                     tempFaceUnknowns,
    double**                                     tempStateSizedVectors,
    double**                                     tempStateSizedSquareMatrices,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  SolverPatch::LimiterStatus limiterStatus = solverPatch.getLimiterStatus();
  if (isRecomputation) {
    limiterStatus = solverPatch.getMergedLimiterStatus(0);

    for (int i=0; i<DIMENSIONS_TIMES_TWO; ++i) {
      assertionEquals(limiterStatus,solverPatch.getMergedLimiterStatus(i));
    } // Dead code elimination will get rid of this loop if Asserts flag is not set.
  }

  logDebug("mergeWithNeighbourDataBasedOnLimiterStatus(...)", "receive data for solver " << _identifier << " from rank " <<
          fromRank << " at vertex x=" << x << ", level=" << level <<
          ", src=" << src << ", dest=" << dest << ",limiterStatus=" << SolverPatch::toString(limiterStatus));

  switch (limiterStatus) {
    case SolverPatch::LimiterStatus::Ok:
    case SolverPatch::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
      _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
      if (!isRecomputation) {
        _solver->mergeWithNeighbourData(
                 fromRank,neighbourTypeAsInt,cellDescriptionsIndex,element,
                 src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      } else {
        _solver->dropNeighbourData(fromRank,src,dest,x,level);
      }
      }break;
    case SolverPatch::LimiterStatus::Troubled:
    case SolverPatch::LimiterStatus::NeighbourIsTroubledCell: {
      const int limiterElement = tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
      assertion(limiterElement!=exahype::solvers::Solver::NotFound);
      _limiter->mergeWithNeighbourData(
          fromRank,neighbourTypeAsInt,cellDescriptionsIndex,limiterElement,
          src,dest,tempFaceUnknowns,tempStateSizedVectors,tempStateSizedSquareMatrices,x,level);
      _solver->dropNeighbourData(fromRank,src,dest,x,level);
      } break;
  }
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMinAndMax(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  const int receivedMinMaxIndex = DataHeap::getInstance().createData(2*_numberOfVariables, 2*_numberOfVariables);
  assertion(DataHeap::getInstance().getData(receivedMinMaxIndex).size()==static_cast<unsigned int>(2*_numberOfVariables));
  double* receivedMinAndMax = DataHeap::getInstance().getData(receivedMinMaxIndex).data();

  DataHeap::getInstance().receiveData(receivedMinAndMax, 2*_numberOfVariables, fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);
  mergeSolutionMinMaxOnFace(solverPatch,faceIndex,receivedMinAndMax,receivedMinAndMax+_numberOfVariables);

  DataHeap::getInstance().deleteData(receivedMinMaxIndex,true);
}

void exahype::solvers::LimitingADERDGSolver::mergeSolutionMinMaxOnFace(
  SolverPatch&  SolverPatch,
  const int           faceIndex,
  const double* const min, const double* const max) const {
  if (SolverPatch.getType() == SolverPatch::Cell ||
      SolverPatch.getType() == SolverPatch::Ancestor ||
      SolverPatch.getType() == SolverPatch::Descendant
      ) {
    double* solutionMin = DataHeap::getInstance().getData( SolverPatch.getSolutionMin()  ).data();
    double* solutionMax = DataHeap::getInstance().getData( SolverPatch.getSolutionMax()  ).data();

    for (int i=0; i<_numberOfVariables; i++) {
      solutionMin[i+faceIndex*_numberOfVariables]  = std::min( solutionMin[i+faceIndex*_numberOfVariables], min[i] );
      solutionMax[i+faceIndex*_numberOfVariables]  = std::max( solutionMax[i+faceIndex*_numberOfVariables], max[i] );
    }
  }
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for(int receives=0; receives<DataMessagesPerNeighbourCommunication; ++receives)
      DataHeap::getInstance().receiveData(
          fromRank, x, level,
          peano::heap::MessageType::NeighbourCommunication);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////
// NEIGHBOUR - Limiter status spreading
///////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendMergedLimiterStatusToNeighbour(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) < dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the "<" !!!
                                                                          // |src|dest| : 1; |dest|src| : 0
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const double mergedLimiterStatusAsDouble =
      static_cast<double>(solverPatch.getMergedLimiterStatus(faceIndex));

  DataHeap::getInstance().sendData(
      &mergedLimiterStatusAsDouble,1,toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataInsteadOfMergedLimiterStatusToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  std::vector<double> emptyMessage(0);
  DataHeap::getInstance().sendData(
      emptyMessage, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithNeighbourMergedLimiterStatus(
    const int                                    fromRank,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
  const int normalOfExchangedFace = tarch::la::equalsReturnIndex(src, dest);
  assertion(normalOfExchangedFace >= 0 && normalOfExchangedFace < DIMENSIONS);
  const int faceIndex = 2 * normalOfExchangedFace +
      (src(normalOfExchangedFace) > dest(normalOfExchangedFace) ? 1 : 0); // !!! Be aware of the ">" !!!

  // Merge the limiter status.
  const int receivedMergedLimiterStatusIndex = DataHeap::getInstance().createData(0, 1);
  assertion(DataHeap::getInstance().getData(receivedMergedLimiterStatusIndex).empty());
  DataHeap::getInstance().receiveData(receivedMergedLimiterStatusIndex,fromRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
  const double* const receivedMergedLimiterStatus =
      DataHeap::getInstance().getData(receivedMergedLimiterStatusIndex).data();
  assertion1(std::lround(*receivedMergedLimiterStatus) <= 3,*receivedMergedLimiterStatus);
  assertion1(std::lround(*receivedMergedLimiterStatus) >= 0,*receivedMergedLimiterStatus);
  SolverPatch::LimiterStatus neighbourMergedLimiterStatus =
      static_cast<SolverPatch::LimiterStatus>(std::lround(*receivedMergedLimiterStatus));
  mergeWithNeighbourLimiterStatus(solverPatch,faceIndex,neighbourMergedLimiterStatus);

  // Clean up.
  DataHeap::getInstance().deleteData(receivedMergedLimiterStatusIndex,true); // TODO(Dominic): We are dealing here with an array of size 1 here.
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourMergedLimiterStatus(
    const int                                    fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) const {
  DataHeap::getInstance().receiveData(fromRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////////////////////////////
// NEIGHBOUR - Solution recomputation
///////////////////////////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendEmptySolverAndLimiterDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  _solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  _limiter->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
}

void exahype::solvers::LimitingADERDGSolver::dropNeighbourSolverAndLimiterData(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) const {
  logDebug("dropNeighbourSolverAndLimiterData(...)", "drop data for solver " << _identifier << " from rank " <<
            fromRank << " at vertex x=" << x << ", level=" << level <<
            ", src=" << src << ", dest=" << dest);

  _limiter->dropNeighbourData(fromRank,src,dest,x,level); // !!! Receive order must be inverted in neighbour comm.
  _solver->dropNeighbourData(fromRank,src,dest,x,level);
}

///////////////////////////////////////
// FORK OR JOIN
///////////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorkerOrMasterDueToForkOrJoin(
      toRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorkerOrMasterDueToForkOrJoin(
          toRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
          toRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
  _limiter->sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
        toRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithWorkerOrMasterDataDueToForkOrJoin(
      fromRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerOrMasterDataDueToForkOrJoin(
        fromRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
        fromRank,x,level);
  } // !!! Receive order must be the same in master<->worker comm.
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerOrMasterDataDueToForkOrJoin(
    const int                                     fromRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
  _limiter->dropWorkerOrMasterDataDueToForkOrJoin(
          fromRank,x,level);
}

///////////////////////////////////
// WORKER->MASTER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                    masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToMaster(masterRank,x,level);
  _limiter->sendDataToMaster(masterRank,x,level);

  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainHasChanged ? 1.0 : -1.0); // TODO(Dominic): ugly
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToMaster(...)","Sending data to master:" <<
             " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(workerRank,x,level); // !!! Receive order must be the same in master<->worker comm.
  _limiter->mergeWithWorkerData(workerRank,x,level);

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool workerLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false;
  updateNextLimiterDomainHasChanged(workerLimiterDomainHasChanged); // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithWorkerData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithWorkerData(...)","_nextLimiterDomainHasChanged=" << _nextLimiterDomainHasChanged);
  }
}

bool exahype::solvers::LimitingADERDGSolver::hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) {
  #if defined(Asserts) || defined(Debug)
  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    assertion(_limiter->hasToSendDataToMaster(cellDescriptionsIndex,limiterElement));
  }
  #endif

  return _solver->hasToSendDataToMaster(cellDescriptionsIndex,element);
}

void exahype::solvers::LimitingADERDGSolver::sendDataToMaster(
    const int                                     masterRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToMaster(
      masterRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToMaster(
        masterRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToMaster(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToMaster(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToMaster(masterRank,x,level);
  _limiter->sendEmptyDataToMaster(masterRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithWorkerData(
    const int                                    workerRank,
    const int                                    workerTypeAsInt,
    const int                                    cellDescriptionsIndex,
    const int                                    element,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithWorkerData(
      workerRank,workerTypeAsInt,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithWorkerData( // !!! Receive order must be the same in master<->worker comm.
        workerRank,workerTypeAsInt,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropWorkerData(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropWorkerData(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropWorkerData(workerRank,x,level);// !!! Receive order must be the same in master<->worker comm.
  _limiter->dropWorkerData(workerRank,x,level);
}

///////////////////////////////////
// MASTER->WORKER
///////////////////////////////////
void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const                                        int workerRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->sendDataToWorker(workerRank,x,level);
  _limiter->sendDataToWorker(workerRank,x,level);

  // TODO(Dominic): Add information that limiter status has been
  // changed for this solver.
  // Send the information to master if limiter status has changed or not
  std::vector<double> dataToSend(0,1);
  dataToSend.push_back(_limiterDomainHasChanged ? 1.0 : -1.0);
  assertion1(dataToSend.size()==1,dataToSend.size());
  if (tarch::parallel::Node::getInstance().getRank()==
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("sendDataToWorker(...)","Sending data to worker:" <<
            " data[0]=" << dataToSend[0]);
  }
  DataHeap::getInstance().sendData(
      dataToSend.data(), dataToSend.size(),
      workerRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const                                        int masterRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  _solver->mergeWithMasterData(masterRank,x,level); // !!! Receive order must be the same in master<->worker comm.
  _limiter->mergeWithMasterData(masterRank,x,level);

  // Receive the information if limiter status has changed
  std::vector<double> receivedData(1); // !!! Creates and fills the vector
  DataHeap::getInstance().receiveData(
      receivedData.data(),receivedData.size(),masterRank, x, level,
      peano::heap::MessageType::MasterWorkerCommunication);
  assertion(tarch::la::equals(receivedData[0],1.0) ||
            tarch::la::equals(receivedData[0],-1.0)); // TODO(Dominic): ugly

  bool masterLimiterDomainHasChanged = tarch::la::equals(receivedData[0],1.0) ? true : false;
  _limiterDomainHasChanged = masterLimiterDomainHasChanged; // !!! It is important that we merge with the "next" field here

  if (tarch::parallel::Node::getInstance().getRank()!=
      tarch::parallel::Node::getInstance().getGlobalMasterRank()) {
    logDebug("mergeWithMasterData(...)","Received data from worker:" <<
            " data[0]=" << receivedData[0]);
    logDebug("mergeWithMasterData(...)","_limiterDomainHasChanged=" << _limiterDomainHasChanged);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendDataToWorker(
    const int                                     workerRank,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendDataToWorker(
      workerRank,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->sendDataToWorker(
        workerRank,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->sendEmptyDataToWorker(workerRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::sendEmptyDataToWorker(
    const int                                     workerRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->sendEmptyDataToWorker(workerRank,x,level);
  _limiter->sendEmptyDataToWorker(workerRank,x,level);
}

void exahype::solvers::LimitingADERDGSolver::mergeWithMasterData(
    const int                                     masterRank,
    const int                                     masterTypeAsInt,
    const int                                     cellDescriptionsIndex,
    const int                                     element,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->mergeWithMasterData(
      masterRank,masterTypeAsInt,cellDescriptionsIndex,element,x,level);

  const int limiterElement = tryGetLimiterElementFromSolverElement(cellDescriptionsIndex,element);
  if (limiterElement!=exahype::solvers::Solver::NotFound) {
    _limiter->mergeWithMasterData( // !!! Receive order must be the same in master<->worker comm.
        masterRank,masterTypeAsInt,cellDescriptionsIndex,limiterElement,x,level);
  } else {
    _limiter->dropMasterData(masterRank,x,level);
  }
}

void exahype::solvers::LimitingADERDGSolver::dropMasterData(
    const int                                     masterRank,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  _solver->dropMasterData(masterRank,x,level);
  _limiter->dropMasterData(masterRank,x,level);
}
#endif

std::string exahype::solvers::LimitingADERDGSolver::toString() const {
  std::ostringstream stringstr;
  toString(stringstr);
  return stringstr.str();
}

void exahype::solvers::LimitingADERDGSolver::toString (std::ostream& out) const {
  out << getIdentifier() << "{_ADERDG: ";
  out << _solver->toString() << "}\n";
  out << getIdentifier() << "{_FV: ";
  out << _limiter->toString() << "}";
}
