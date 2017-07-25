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

#include "exahype/mappings/InitialCondition.h"

#include <cmath>

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "tarch/multicore/Loop.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/TimeStepSizeComputation.h"

void exahype::mappings::InitialCondition::prepareTemporaryVariables() {
  assertion(_limiterDomainHasChanged ==nullptr);
  int numberOfSolvers      = exahype::solvers::RegisteredSolvers.size();
  _limiterDomainHasChanged = new bool    [numberOfSolvers];

  for (unsigned int solverNumber=0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    _limiterDomainHasChanged[solverNumber] = false;
  }
}

void exahype::mappings::InitialCondition::deleteTemporaryVariables() {
  if (_limiterDomainHasChanged!=nullptr) {
    delete[] _limiterDomainHasChanged;
    _limiterDomainHasChanged = nullptr;
  }
}

peano::CommunicationSpecification
exahype::mappings::InitialCondition::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
      MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
      MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::InitialCondition::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::InitialCondition::_log(
    "exahype::mappings::InitialCondition");

void exahype::mappings::InitialCondition::enterCell(
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
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined1);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
      auto solver = exahype::solvers::RegisteredSolvers[i];

      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->prepareNextNeighbourMerging(
            fineGridCell.getCellDescriptionsIndex(),element,
            fineGridVertices,fineGridVerticesEnumerator);

        solver->setInitialConditions(
            fineGridCell.getCellDescriptionsIndex(),
            element,
            fineGridVertices,
            fineGridVerticesEnumerator);

        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
          bool limiterDomainHasChanged =
              static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->
              updateMergedLimiterStatusAndMinAndMaxAfterSetInitialConditions(fineGridCell.getCellDescriptionsIndex(),element);
          _limiterDomainHasChanged[i] |= limiterDomainHasChanged;
        }
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::InitialCondition::InitialCondition(
    const InitialCondition& masterThread)
    :
    _limiterDomainHasChanged(nullptr) {
  prepareTemporaryVariables();
}
void exahype::mappings::InitialCondition::mergeWithWorkerThread(
    const InitialCondition& workerThread) {
  for (int i = 0; i < static_cast<int>(exahype::solvers::RegisteredSolvers.size()); i++) {
    _limiterDomainHasChanged[i] |= workerThread._limiterDomainHasChanged[i];
  }
}
#endif

void exahype::mappings::InitialCondition::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  prepareTemporaryVariables();

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::InitialCondition::endIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("endIteration(State)", solverState);

  // 1. Merge limiter domain status
  bool limiterDomainHasChanged = solverState.limiterDomainHasChanged();
  for (unsigned int solverNumber = 0; solverNumber < exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    if (exahype::solvers::RegisteredSolvers[solverNumber]->getType()==exahype::solvers::Solver::Type::LimitingADERDG) {
      logInfo("endIteration(State)", "_limiterDomainHasChanged[solverNumber]="<<_limiterDomainHasChanged[solverNumber]<<
                    ",limiterDomainHasChanged="<<limiterDomainHasChanged);

      auto* limitingADERDGSolver =
          static_cast<exahype::solvers::LimitingADERDGSolver*>(exahype::solvers::RegisteredSolvers[solverNumber]);
      limitingADERDGSolver->updateNextLimiterDomainHasChanged(_limiterDomainHasChanged[solverNumber]);

      limiterDomainHasChanged |= limitingADERDGSolver->getNextLimiterDomainHasChanged();

      logInfo("endIteration(State)", "limitingADERDGSolver->getNextLimiterDomainHasChanged()="<<limitingADERDGSolver->getNextLimiterDomainHasChanged()<<",_limiterDomainHasChanged[solverNumber]="<<_limiterDomainHasChanged[solverNumber]<<
              ",limiterDomainHasChanged="<<limiterDomainHasChanged);
    }
  }
  solverState.updateLimiterDomainHasChanged(limiterDomainHasChanged);

  // 2.
  deleteTemporaryVariables();
  // 3. Veto the fused time stepping time step size reinitialisation in TimeStepSizeComputation::enterCell
  exahype::mappings::TimeStepSizeComputation::VetoFusedTimeSteppingTimeStepSizeReinitialisation = true;

  logTraceOutWith1Argument("endIteration(State)", solverState);
}

//
// Below all methods are nop.
//
//===================================

#ifdef Parallel
void exahype::mappings::InitialCondition::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::InitialCondition::prepareSendToWorker(
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

void exahype::mappings::InitialCondition::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialCondition::mergeWithMaster(
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

void exahype::mappings::InitialCondition::receiveDataFromMaster(
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

void exahype::mappings::InitialCondition::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::InitialCondition::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::InitialCondition::InitialCondition()
    :
    _limiterDomainHasChanged(nullptr){
  // do nothing
}

exahype::mappings::InitialCondition::~InitialCondition() {
  // do nothing
}

void exahype::mappings::InitialCondition::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialCondition::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::InitialCondition::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::InitialCondition::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::InitialCondition::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::InitialCondition::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
