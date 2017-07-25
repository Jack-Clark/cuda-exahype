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

#include "exahype/mappings/Sending.h"

#include <limits>

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/Lock.h"
#include "tarch/parallel/Node.h"

#include "peano/utils/Globals.h"
#include "peano/utils/Loop.h"

#include "peano/datatraversal/autotuning/Oracle.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/VertexOperations.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"



tarch::multicore::BooleanSemaphore exahype::mappings::Sending::_semaphoreForRestriction;

#ifdef Parallel
bool exahype::mappings::Sending::SkipReductionInBatchedTimeSteps = false;
#endif

peano::CommunicationSpecification
exahype::mappings::Sending::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::SendDataAndStateBeforeFirstTouchVertexFirstTime,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::SendDataAndStateAfterProcessingOfLocalSubtree,
      true);
}

/**
 * Nop.
 */
peano::MappingSpecification exahype::mappings::Sending::
    touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification exahype::mappings::Sending::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Sending::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Sending::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Sending::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::Sending::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::Sending::_log(
    "exahype::mappings::Sending");

void exahype::mappings::Sending::beginIteration(
    exahype::State& solverState) {
  _localState = solverState;

  #ifdef Parallel
  if (_localState.getSendMode()!=exahype::records::State::SendMode::SendNothing) {
    exahype::solvers::ADERDGSolver::Heap::getInstance().startToSendSynchronousData();
    exahype::solvers::FiniteVolumesSolver::Heap::getInstance().startToSendSynchronousData();
    DataHeap::getInstance().startToSendSynchronousData();
    MetadataHeap::getInstance().startToSendSynchronousData();
  }
  #endif

  logDebug("beginIteration(...)","MergeMode="<<_localState.getMergeMode()<<", SendMode="<<_localState.getSendMode());
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Sending::Sending(const Sending& masterThread) :
  _localState(masterThread._localState) {}

void exahype::mappings::Sending::mergeWithWorkerThread(
    const Sending& workerThread) {
  // do nothing for now
}
#endif

void exahype::mappings::Sending::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (_localState.getSendMode()==exahype::records::State::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::ReduceAndMergeTimeStepDataAndSendFaceData) {
    const int numberOfSolvers = static_cast<int>(exahype::solvers::RegisteredSolvers.size());

    auto grainSize = peano::datatraversal::autotuning::Oracle::
        getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined12);
    pfor(solverNumber, 0, numberOfSolvers, grainSize.getGrainSize())
      exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        solver->prolongateDataAndPrepareDataRestriction(fineGridCell.getCellDescriptionsIndex(),element);
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }
}


void exahype::mappings::Sending::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("leaveCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (_localState.getSendMode()==exahype::records::State::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::ReduceAndMergeTimeStepDataAndSendFaceData) {
    const int numberOfSolvers = static_cast<int>(exahype::solvers::RegisteredSolvers.size());
    auto grainSize = peano::datatraversal::autotuning::Oracle::
        getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined13);
    pfor(solverNumber, 0, numberOfSolvers, grainSize.getGrainSize())
      exahype::solvers::Solver* solver = exahype::solvers::RegisteredSolvers[solverNumber];

      int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);
      if (element!=exahype::solvers::Solver::NotFound) {
        exahype::solvers::Solver* solver =
            exahype::solvers::RegisteredSolvers[solverNumber];

        exahype::solvers::Solver::SubcellPosition subcellPosition =
            solver->computeSubcellPositionOfCellOrAncestor(
                fineGridCell.getCellDescriptionsIndex(),element);

        if (subcellPosition.parentElement!=exahype::solvers::Solver::NotFound &&
            exahype::amr::onBoundaryOfParent(
                subcellPosition.subcellIndex,subcellPosition.levelDifference)) {
          tarch::multicore::Lock lock(_semaphoreForRestriction);

          solver->restrictData(fineGridCell.getCellDescriptionsIndex(),
                               element,
                               subcellPosition.parentCellDescriptionsIndex,
                               subcellPosition.parentElement,
                               subcellPosition.subcellIndex);
        }

        solver->postProcess(fineGridCell.getCellDescriptionsIndex(),element);
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }

  logTraceOutWith1Argument("leaveCell(...)", fineGridCell);
}


#ifdef Parallel
///////////////////////////////////////
// NEIGHBOUR
///////////////////////////////////////
void exahype::mappings::Sending::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  if (_localState.getSendMode()==exahype::records::State::SendMode::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData) {
    dfor2(dest)
      dfor2(src)
      if (vertex.hasToSendMetadata(src,dest,toRank)) {
        vertex.tryDecrementFaceDataExchangeCountersOfSource(src,dest);
        if (vertex.hasToSendDataToNeighbour(src,dest)) {
          sendSolverDataToNeighbour(
              toRank,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              x,level);
        } else {
          sendEmptySolverDataToNeighbour(toRank,src,dest,x,level);
        }
      }
      enddforx
    enddforx
  }
}

void exahype::mappings::Sending::sendEmptySolverDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS, int>&     src,
    const tarch::la::Vector<DIMENSIONS, int>&     dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level) {
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
  }

  auto encodedMetadata = exahype::Vertex::createEncodedMetadataSequenceWithInvalidEntries();
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

void exahype::mappings::Sending::sendSolverDataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level) {
  assertion(exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));
  assertion(exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(srcCellDescriptionIndex));

  int solverNumber=0;
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    int element = solver->tryGetElement(srcCellDescriptionIndex,solverNumber);
    if (element!=exahype::solvers::Solver::NotFound) {
      solver->sendDataToNeighbour(toRank,srcCellDescriptionIndex,element,src,dest,x,level);
    } else {
      solver->sendEmptyDataToNeighbour(toRank,src,dest,x,level);
    }
    ++solverNumber;
  }

  auto encodedMetadata = exahype::Vertex::encodeMetadata(srcCellDescriptionIndex);
  MetadataHeap::getInstance().sendData(
      encodedMetadata, toRank, x, level,
      peano::heap::MessageType::NeighbourCommunication);
}

///////////////////////////////////////
// WORKER->MASTER->WORKER
///////////////////////////////////////
void exahype::mappings::Sending::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  if (_localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData) {
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      solver->sendDataToMaster(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());
    }
  }

  if (_localState.getSendMode()==exahype::records::State::SendMode::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData) {
    if (localCell.isInside() && localCell.isInitialised()) {
      exahype::Vertex::sendEncodedMetadata(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          localCell.getCellDescriptionsIndex(),
          peano::heap::MessageType::MasterWorkerCommunication,
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());

      int solverNumber=0;
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        int element = solver->tryGetElement(localCell.getCellDescriptionsIndex(),solverNumber);
        if (element!=exahype::solvers::Solver::NotFound) {
          solver->sendDataToMaster(
              tarch::parallel::NodePool::getInstance().getMasterRank(),
              localCell.getCellDescriptionsIndex(),
              element,
              verticesEnumerator.getCellCenter(),
              verticesEnumerator.getLevel());
        } else {
          solver->sendEmptyDataToMaster(
              tarch::parallel::NodePool::getInstance().getMasterRank(),
              verticesEnumerator.getCellCenter(),
              verticesEnumerator.getLevel());
        }
        ++solverNumber;
      }
    } else if (localCell.isInside() && !localCell.isInitialised()) {
      exahype::Vertex::sendEncodedMetadataSequenceWithInvalidEntries(
          tarch::parallel::NodePool::getInstance().getMasterRank(),
          peano::heap::MessageType::MasterWorkerCommunication,
          verticesEnumerator.getCellCenter(),
          verticesEnumerator.getLevel());

      for (auto solver : exahype::solvers::RegisteredSolvers) {
        solver->sendEmptyDataToMaster(
            tarch::parallel::NodePool::getInstance().getMasterRank(),
            verticesEnumerator.getCellCenter(),
            verticesEnumerator.getLevel());
      }
    } // else  do nothing
  }
}

void exahype::mappings::Sending::mergeWithMaster(
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
  if (_localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData) {
    for (auto& p : exahype::solvers::RegisteredSolvers) {
      p->mergeWithWorkerData(
          worker,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());
    }
  }

  if (_localState.getSendMode()==exahype::records::State::SendMode::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData) {
    // TODO(Dominic): Add to docu, we only know from the worker grid cell if it is inside.
    // I encountered a problem where workerGridCell.isInside() != fineGridCell.isInside()
    if (workerGridCell.isInside() && fineGridCell.isInitialised()) {
      int receivedMetadataIndex = exahype::MetadataHeap::getInstance().createData(
          0,exahype::solvers::RegisteredSolvers.size());
      exahype::MetadataHeap::getInstance().receiveData(
          receivedMetadataIndex,worker,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel(),
          peano::heap::MessageType::MasterWorkerCommunication);
      MetadataHeap::HeapEntries& receivedMetadata =
          MetadataHeap::getInstance().getData(receivedMetadataIndex);

      int solverNumber=0;
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

        if (element!=exahype::solvers::Solver::NotFound &&
            receivedMetadata[solverNumber].getU()!=exahype::Vertex::InvalidMetadataEntry) {
          solver->mergeWithWorkerData(
              worker,
              receivedMetadata[solverNumber].getU(),
              fineGridCell.getCellDescriptionsIndex(),element,
              fineGridVerticesEnumerator.getCellCenter(),
              fineGridVerticesEnumerator.getLevel());
        } else {
          solver->dropWorkerData(
              worker,
              fineGridVerticesEnumerator.getCellCenter(),
              fineGridVerticesEnumerator.getLevel());
        }

        ++solverNumber;
      }
      exahype::MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
    } else if (workerGridCell.isInside() && !fineGridCell.isInitialised()) {
      exahype::Vertex::dropMetadata(
          worker,
          peano::heap::MessageType::MasterWorkerCommunication,
          fineGridVerticesEnumerator.getCellCenter(),
          fineGridVerticesEnumerator.getLevel());

      for (auto solver : exahype::solvers::RegisteredSolvers) {
        solver->dropWorkerData(
            worker,
            fineGridVerticesEnumerator.getCellCenter(),
            fineGridVerticesEnumerator.getLevel());
      }
    } // else do nothing
  }
}


bool exahype::mappings::Sending::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  bool workerHasToSendDataMaster = false;

  if ((_localState.getSendMode()==exahype::records::State::SendMode::SendFaceData ||
      _localState.getSendMode()==exahype::records::State::SendMode::ReduceAndMergeTimeStepDataAndSendFaceData)
      && fineGridCell.isInside()) {
    if (fineGridCell.isInitialised()) {
      int solverNumber=0;
      for (auto solver : exahype::solvers::RegisteredSolvers) {
        int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),solverNumber);

        if (element!=exahype::solvers::Solver::NotFound) {
          workerHasToSendDataMaster |=
              solver->hasToSendDataToMaster(fineGridCell.getCellDescriptionsIndex(),element);
        }

        ++solverNumber;
      }
    }
  }

  if (!peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated()
      &&
      workerHasToSendDataMaster
      &&
      SkipReductionInBatchedTimeSteps) {
    return false;
  }
  else return true;
}

//
// Below all functions are nop.
//
// ====================================

void exahype::mappings::Sending::receiveDataFromMaster(
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

void exahype::mappings::Sending::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Sending::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Sending::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Sending::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Sending::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Sending::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Sending::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}
#endif

exahype::mappings::Sending::Sending() {
  // do nothing
}

exahype::mappings::Sending::~Sending() {
  // do nothing
}

void exahype::mappings::Sending::endIteration(
    exahype::State& solverState) {
  // do nothing
}

void exahype::mappings::Sending::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Sending::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Sending::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Sending::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Sending::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
