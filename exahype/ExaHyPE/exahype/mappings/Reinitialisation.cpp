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

#include "exahype/mappings/Reinitialisation.h"

#include <cmath>

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "tarch/multicore/Loop.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/LimiterStatusSpreading.h"

peano::CommunicationSpecification
exahype::mappings::Reinitialisation::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
      MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
      MaskOutWorkerMasterDataAndStateExchange,
      true);
}

// Everything below is nop.
peano::MappingSpecification
exahype::mappings::Reinitialisation::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::Reinitialisation::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial,true);
}

peano::MappingSpecification
exahype::mappings::Reinitialisation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::Reinitialisation::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::Reinitialisation::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::Reinitialisation::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::Reinitialisation::_log(
    "exahype::mappings::Reinitialisation");

void exahype::mappings::Reinitialisation::beginIteration(
    exahype::State& solverState) {

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceMerges = 0;
  _boundaryFaceMerges = 0;
  #endif
}

void exahype::mappings::Reinitialisation::endIteration(
    exahype::State& solverState) {
  // TODO(Dominic): Assess

  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      limitingADERDGSolver->rollbackToPreviousTimeStep();

      if (!exahype::State::fuseADERDGPhases()) {
        limitingADERDGSolver->reconstructStandardTimeSteppingDataAfterRollback();
      }
    }
  }

  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interior face merges: " << _interiorFaceMerges);
  logInfo("endIteration(...)","boundary face merges: " << _boundaryFaceMerges);
  #endif
}

void exahype::mappings::Reinitialisation::enterCell(
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
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined10);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
      auto solver = exahype::solvers::RegisteredSolvers[i];

      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
      if (element!=exahype::solvers::Solver::NotFound) {
        if(solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
           && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
          auto limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

          limitingADERDGSolver->updateMergedLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element); // update before reinitialisation

          limitingADERDGSolver->rollbackToPreviousTimeStep(fineGridCell.getCellDescriptionsIndex(),element); // Loads back the old corrector time step size.
          if (!exahype::State::fuseADERDGPhases()) {
            limitingADERDGSolver->reconstructStandardTimeSteppingDataAfterRollback(
                fineGridCell.getCellDescriptionsIndex(),element);
          }

          limitingADERDGSolver->reinitialiseSolvers(fineGridCell.getCellDescriptionsIndex(),element,
              fineGridCell,fineGridVertices,fineGridVerticesEnumerator); // TODO(Dominic): Probably need to merge those
        }

        solver->prepareNextNeighbourMerging(
            fineGridCell.getCellDescriptionsIndex(),element,
            fineGridVertices,fineGridVerticesEnumerator); // !!! Has to be done after reinitialisation since we might add new finite volumes patches here.
                                                          // !!! Has to be done for all solvers (cf. touchVertexFirstTime etc.)
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::Reinitialisation::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  logTraceInWith6Arguments("touchVertexFirstTime(...)", fineGridVertex,
                             fineGridX, fineGridH,
                             coarseGridVerticesEnumerator.toString(),
                             coarseGridCell, fineGridPositionOfVertex);

  dfor2(pos1)
    dfor2(pos2)
      if (fineGridVertex.hasToMergeNeighbours(pos1,pos2)) { // Assumes that we have to valid indices // TODO(Dominic): Probably have to consider Voronoi neighbours later on when we use high order schemes
        auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined11);
        pfor(solverNumber, 0, static_cast<int>(solvers::RegisteredSolvers.size()),grainSize.getGrainSize())
          auto solver = exahype::solvers::RegisteredSolvers[solverNumber];

          if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
              && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
            const int cellDescriptionsIndex1 = fineGridVertex.getCellDescriptionsIndex()[pos1Scalar];
            const int cellDescriptionsIndex2 = fineGridVertex.getCellDescriptionsIndex()[pos2Scalar];
            const int element1 = solver->tryGetElement(cellDescriptionsIndex1,solverNumber);
            const int element2 = solver->tryGetElement(cellDescriptionsIndex2,solverNumber);
            if (element2>=0 && element1>=0) {
              auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
              limitingADERDGSolver->mergeLimiterStatusOfNeighbours(
                  cellDescriptionsIndex1,element1,cellDescriptionsIndex2,element2,pos1,pos2);
            }
          }

          #ifdef Debug // TODO(Dominic)
          _interiorFaceMerges++;
          #endif
        endpfor
        grainSize.parallelSectionHasTerminated();

        fineGridVertex.setMergePerformed(pos1,pos2,true);
      }
    enddforx
  enddforx

  logTraceOutWith1Argument("touchVertexFirstTime(...)", fineGridVertex);
}

#ifdef Parallel
void exahype::mappings::Reinitialisation::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  dfor2(dest)
    dfor2(src)
      if (vertex.hasToSendMetadata(src,dest,toRank)) {
        vertex.tryDecrementFaceDataExchangeCountersOfSource(src,dest);
        if (vertex.hasToSendDataToNeighbour(src,dest)) { // Only comm. data once per face
          sendDataToNeighbour(
              toRank,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              x,level);
        } else {
          sendEmptyDataToNeighbour(
              toRank,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              x,level);
        }
      }
    enddforx
  enddforx
}

/*
 * We only send empty data for LimitingADERDGSolvers
 * where we have detected a change of the limiter domain.
 * This information should be available on all ranks.
 * We ignore other solver types.
 */
void exahype::mappings::Reinitialisation::sendEmptyDataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS, int>&    src,
    const tarch::la::Vector<DIMENSIONS, int>&    dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  for (unsigned   int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {

      logDebug("sendEmptyDataToNeighbour(...)", "send empty data for solver " << solverNumber << " to rank " <<
              toRank << " at vertex x=" << x << ", level=" << level <<
              ", src=" << src << ", dest=" << dest);

      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      limitingADERDGSolver->sendEmptySolverAndLimiterDataToNeighbour(toRank,src,dest,x,level);
    }
  }

  auto encodedMetadata = exahype::Vertex::createEncodedMetadataSequenceWithInvalidEntries();
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

/*
 * We only send empty data for LimitingADERDGSolvers
 * where we have detected a change of the limiter domain.
 * This information should be available on all ranks.
 * We ignore other solver types.
 */
void exahype::mappings::Reinitialisation::sendDataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));

  for (unsigned   int solverNumber=0; solverNumber<exahype::solvers::RegisteredSolvers.size(); ++solverNumber) {
    auto* solver = exahype::solvers::RegisteredSolvers[solverNumber];

    if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
        && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
      auto* limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);
      const int element = solver->tryGetElement(srcCellDescriptionIndex,solverNumber);

      if (element!=exahype::solvers::Solver::NotFound) { // Check if a patch exists on the cell

        logDebug("sendDataToNeighbour(...)", "send data for solver " << solverNumber << " to rank " <<
                toRank << " at vertex x=" << x << ", level=" << level <<
                ", src=" << src << ", dest=" << dest);

        limitingADERDGSolver->sendDataToNeighbourBasedOnLimiterStatus(
            toRank,srcCellDescriptionIndex,element,src,dest,
            true, /* isRecomputation */
            x,level);
      } else {

        logDebug("sendDataToNeighbour(...)", "send empty data for solver " << solverNumber << " to rank " <<
                 toRank << " at vertex x=" << x << ", level=" << level <<
                 ", src=" << src << ", dest=" << dest);

        limitingADERDGSolver->sendEmptySolverAndLimiterDataToNeighbour(
            toRank,src,dest,x,level);
      }
    }
  }

  auto encodedMetadata = exahype::Vertex::encodeMetadata(srcCellDescriptionIndex);
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}



//
// Below all methods are nop.
//
//===================================



void exahype::mappings::Reinitialisation::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Reinitialisation::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Reinitialisation::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Reinitialisation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Reinitialisation::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::Reinitialisation::prepareSendToWorker(
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

void exahype::mappings::Reinitialisation::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Reinitialisation::mergeWithMaster(
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

void exahype::mappings::Reinitialisation::receiveDataFromMaster(
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

void exahype::mappings::Reinitialisation::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Reinitialisation::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::Reinitialisation::Reinitialisation()
  #ifdef Debug
  :
  _interiorFaceMerges(0),
  _boundaryFaceMerges(0)
  #endif
{
  // do nothing
}

exahype::mappings::Reinitialisation::~Reinitialisation() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Reinitialisation::Reinitialisation(
    const Reinitialisation& masterThread)
#ifdef Debug
  :
  _interiorFaceMerges(0),
  _boundaryFaceMerges(0)
  #endif
{
  // do nothing
}
void exahype::mappings::Reinitialisation::mergeWithWorkerThread(
    const Reinitialisation& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::Reinitialisation::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Reinitialisation::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Reinitialisation::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Reinitialisation::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Reinitialisation::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Reinitialisation::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
