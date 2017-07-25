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
 
#include "exahype/mappings/SolutionUpdate.h"

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "tarch/multicore/Loop.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/VertexOperations.h"
#include "multiscalelinkedcell/HangingVertexBookkeeper.h"


peano::CommunicationSpecification
exahype::mappings::SolutionUpdate::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

/**
 * @todo Please tailor the parameters to your mapping's properties.
 */
peano::MappingSpecification
exahype::mappings::SolutionUpdate::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}


peano::MappingSpecification
exahype::mappings::SolutionUpdate::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::SolutionUpdate::_log(
    "exahype::mappings::SolutionUpdate");

void exahype::mappings::SolutionUpdate::prepareTemporaryVariables() {
  assertion(_limiterDomainHasChanged ==nullptr);
  assertion(_tempUnknowns            ==nullptr);
  assertion(_tempStateSizedVectors   ==nullptr);

  int numberOfSolvers      = exahype::solvers::RegisteredSolvers.size();
  _limiterDomainHasChanged = new bool    [numberOfSolvers];
  _tempStateSizedVectors   = new double**[numberOfSolvers];
  _tempUnknowns            = new double**[numberOfSolvers];

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    if  (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes ||
        solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      const int numberOfStateSizedVectors = 1+2*DIMENSIONS; // max; see riemannSolverNonlinear(5) or kernels::finitevolumes::godunov::solutionUpdate (1+2*DIMENSIONS)
      _tempStateSizedVectors[solverNumber]    = new double*[numberOfStateSizedVectors];
      _tempStateSizedVectors[solverNumber][0] = new double[numberOfStateSizedVectors * solver->getNumberOfVariables()];
      for (int i=1; i<numberOfStateSizedVectors; ++i) {
        _tempStateSizedVectors[solverNumber][i] = _tempStateSizedVectors[solverNumber][i-1]+solver->getNumberOfVariables();
      }
      //
      // TODO(Dominic): This will change if we use a different method than a 1st order Godunov method:
      _tempUnknowns[solverNumber] = nullptr;
    } else {
      _tempStateSizedVectors[solverNumber] = nullptr;
      _tempUnknowns     [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::mappings::SolutionUpdate::prepareLimiterDomainHasChangedFlags() {
  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _limiterDomainHasChanged[solverNumber] = false;
  }
}

void exahype::mappings::SolutionUpdate::deleteTemporaryVariables() {
  if (_tempStateSizedVectors!=nullptr) {
    assertion(_limiterDomainHasChanged!=nullptr);
    assertion(_tempUnknowns           !=nullptr);

    int solverNumber=0;
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      if (solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes ||
          solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
        //
        delete[] _tempStateSizedVectors[solverNumber][0];
        delete[] _tempStateSizedVectors[solverNumber];
        _tempStateSizedVectors[solverNumber] = nullptr;
        // TODO(Dominic): This will change if we use a different method than a 1st order Godunov method:
        _tempUnknowns[solverNumber] = nullptr;
      }

      ++solverNumber;
    }

    delete[] _tempStateSizedVectors;
    delete[] _tempUnknowns;
    delete[] _limiterDomainHasChanged;
    _tempStateSizedVectors   = nullptr;
    _tempUnknowns            = nullptr;
    _limiterDomainHasChanged = nullptr;
  }
}

exahype::mappings::SolutionUpdate::SolutionUpdate() {
  // do nothing
}

exahype::mappings::SolutionUpdate::~SolutionUpdate() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::SolutionUpdate::SolutionUpdate(
    const SolutionUpdate& masterThread)
: _tempStateSizedVectors(nullptr),
  _tempUnknowns(nullptr),
  _limiterDomainHasChanged(nullptr) {
  prepareTemporaryVariables();
  prepareLimiterDomainHasChangedFlags();
}

void exahype::mappings::SolutionUpdate::mergeWithWorkerThread(
    const SolutionUpdate& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _limiterDomainHasChanged[i] |= workerThread._limiterDomainHasChanged[i];
  }
}
#endif

void exahype::mappings::SolutionUpdate::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfSolvers = exahype::solvers::RegisteredSolvers.size();
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined17);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
      auto solver = exahype::solvers::RegisteredSolvers[i];

      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->updateSolution(
            fineGridCell.getCellDescriptionsIndex(),
            element,
            _tempStateSizedVectors[i],
            _tempUnknowns[i],
            fineGridVertices,
            fineGridVerticesEnumerator);

        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          bool limiterDomainHasChanged =
              static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
              updateMergedLimiterStatusAndMinAndMaxAfterSolutionUpdate(fineGridCell.getCellDescriptionsIndex(),element);
          _limiterDomainHasChanged[i] |= limiterDomainHasChanged;
        }

        solver->prepareNextNeighbourMerging(
            fineGridCell.getCellDescriptionsIndex(),element,
            fineGridVertices,fineGridVerticesEnumerator);
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::SolutionUpdate::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  prepareTemporaryVariables();
  prepareLimiterDomainHasChangedFlags();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::SolutionUpdate::endIteration(
    exahype::State& solverState) { // TODO(Dominic): Can't stay here
  logTraceInWith1Argument("endIteration(State)", solverState);

  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    if (exahype::solvers::RegisteredSolvers[solverNumber]->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      auto* limitingADERDGSolver =
          static_cast<exahype::solvers::LimitingADERDGSolver*>(exahype::solvers::RegisteredSolvers[solverNumber]);
      limitingADERDGSolver->updateNextLimiterDomainHasChanged(_limiterDomainHasChanged[solverNumber]);

      logDebug("endIteration(State)", "solver "<<solverNumber<<": next limiter domain has changed: "<<limitingADERDGSolver->getNextLimiterDomainHasChanged());
    }
  }

  deleteTemporaryVariables();

  logTraceOutWith1Argument("endIteration(State)", solverState);
}



//
// Below all methods are nop.
//
//=====================================



void exahype::mappings::SolutionUpdate::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

#ifdef Parallel
void exahype::mappings::SolutionUpdate::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::SolutionUpdate::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}

void exahype::mappings::SolutionUpdate::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithMaster(
    const exahype::Cell& workerGridCell,
    exahype::Vertex* const workerGridVertices,
    const peano::grid::VertexEnumerator& workerEnumerator,
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker, const exahype::State& workerState,
    exahype::State& masterState) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::SolutionUpdate::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::SolutionUpdate::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
