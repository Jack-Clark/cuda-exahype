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

#include "exahype/mappings/LimiterStatusMergingAndSpreadingMPI.h"

#include <cmath>

#include "peano/utils/Globals.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "tarch/multicore/Loop.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/mappings/LimiterStatusSpreading.h"

#include "exahype/mappings/LimiterStatusMergingMPI.h"

peano::CommunicationSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
      MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
      MaskOutWorkerMasterDataAndStateExchange,
      true);
}

// Everything below is nop.
peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::Serial,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::LimiterStatusMergingAndSpreadingMPI::_log(
    "exahype::mappings::LimiterStatusMergingAndSpreadingMPI");

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::enterCell(
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
    auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(numberOfSolvers, peano::datatraversal::autotuning::MethodTrace::UserDefined2);
    pfor(i, 0, numberOfSolvers, grainSize.getGrainSize())
      auto solver = exahype::solvers::RegisteredSolvers[i];

      const int element = solver->tryGetElement(fineGridCell.getCellDescriptionsIndex(),i);
      if (element!=exahype::solvers::Solver::NotFound) {
        if (solver->getType()==exahype::solvers::Solver::Type::LimitingADERDG
            && static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiterDomainHasChanged()) {
          auto limitingADERDGSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver);

          limitingADERDGSolver->updateMergedLimiterStatus(fineGridCell.getCellDescriptionsIndex(),element);
        }

        solver->prepareNextNeighbourMerging(
            fineGridCell.getCellDescriptionsIndex(),element,
            fineGridVertices,fineGridVerticesEnumerator); // !!! Has to be done for all solvers (cf. touchVertexFirstTime etc.)
      }
    endpfor
    grainSize.parallelSectionHasTerminated();
  }

  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexFirstTime(
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
            parallelise(solvers::RegisteredSolvers.size(), peano::datatraversal::autotuning::MethodTrace::UserDefined3);
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

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::beginIteration(
    exahype::State& solverState) {
  #ifdef Parallel
  DataHeap::getInstance().startToSendSynchronousData();
  #endif

  #ifdef Debug // TODO(Dominic): And not parallel and not shared memory
  _interiorFaceMerges = 0;
  _boundaryFaceMerges = 0;
  #endif
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::endIteration(
    exahype::State& solverState) {
  #if defined(Debug) // TODO(Dominic): Use logDebug if it works with filters
  logInfo("endIteration(...)","interior face merges: " << _interiorFaceMerges); /// @deprecated fields
  logInfo("endIteration(...)","boundary face merges: " << _boundaryFaceMerges);
  #endif
}

#ifdef Parallel
void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  dfor2(myDest)
    dfor2(mySrc)
      tarch::la::Vector<DIMENSIONS, int> dest = tarch::la::Vector<DIMENSIONS, int>(1) - myDest;
      tarch::la::Vector<DIMENSIONS, int> src  = tarch::la::Vector<DIMENSIONS, int>(1) - mySrc;

      int destScalar = TWO_POWER_D - myDestScalar - 1;
      int srcScalar  = TWO_POWER_D - mySrcScalar  - 1;

      if (vertex.hasToReceiveMetadata(src,dest,fromRank)) {
        int receivedMetadataIndex = MetadataHeap::getInstance().
            createData(0,exahype::solvers::RegisteredSolvers.size());
        MetadataHeap::getInstance().receiveData(
            receivedMetadataIndex,
            fromRank, fineGridX, level,
            peano::heap::MessageType::NeighbourCommunication);
        exahype::MetadataHeap::HeapEntries& receivedMetadata = MetadataHeap::getInstance().getData(receivedMetadataIndex);
        assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());

        if(vertex.hasToMergeWithNeighbourData(src,dest)) { // Only comm. data once per face
          exahype::mappings::LimiterStatusMergingMPI::mergeNeighourMergedLimiterStatus(
              fromRank,
              src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              fineGridX,level,
              receivedMetadata);

          vertex.setFaceDataExchangeCountersOfDestination(src,dest,TWO_POWER_D); // !!! Do not forget this
          vertex.setMergePerformed(src,dest,true);
        } else {
          exahype::mappings::LimiterStatusMergingMPI::dropNeighbourMergedLimiterStatus(
              fromRank,
              src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              fineGridX,level,
              receivedMetadata);
        }
        // Clean up
        MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
      }
    enddforx
  enddforx
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  dfor2(dest)
    dfor2(src)
      if (vertex.hasToSendMetadata(src,dest,toRank)) {
        vertex.tryDecrementFaceDataExchangeCountersOfSource(src,dest);
        if (vertex.hasToSendDataToNeighbour(src,dest)) {
          exahype::mappings::LimiterStatusSpreading::sendMergedLimiterStatusToNeighbour(
              toRank,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              x,level);
        } else {
          exahype::mappings::LimiterStatusSpreading::sendEmptyDataInsteadOfMergedLimiterStatusToNeighbour(
              toRank,src,dest,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              x,level);
        }
      }
    enddforx
  enddforx
}


//
// Below all methods are nop.
//
//===================================

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

bool exahype::mappings::LimiterStatusMergingAndSpreadingMPI::prepareSendToWorker(
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

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithMaster(
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

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::receiveDataFromMaster(
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

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

exahype::mappings::LimiterStatusMergingAndSpreadingMPI::LimiterStatusMergingAndSpreadingMPI()
#ifdef Debug
:
_interiorFaceMerges(0),
_boundaryFaceMerges(0)
#endif
{
  // do nothing
}

exahype::mappings::LimiterStatusMergingAndSpreadingMPI::~LimiterStatusMergingAndSpreadingMPI() {
  // do nothing
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::LimiterStatusMergingAndSpreadingMPI::LimiterStatusMergingAndSpreadingMPI(
    const LimiterStatusMergingAndSpreadingMPI& masterThread)
  #ifdef Debug
  :
  _interiorFaceMerges(0),
  _boundaryFaceMerges(0)
  #endif
{
  // do nothing
}
void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::mergeWithWorkerThread(
    const LimiterStatusMergingAndSpreadingMPI& workerThread) {
  // do nothing
}
#endif

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::LimiterStatusMergingAndSpreadingMPI::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
