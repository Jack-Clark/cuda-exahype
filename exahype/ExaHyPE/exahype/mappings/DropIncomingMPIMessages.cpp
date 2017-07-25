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
 
#include "exahype/mappings/DropIncomingMPIMessages.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/Solver.h"

peano::CommunicationSpecification
exahype::mappings::DropIncomingMPIMessages::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::
          MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::
          MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::DropIncomingMPIMessages::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification exahype::mappings::DropIncomingMPIMessages::
    touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::DropIncomingMPIMessages::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::DropIncomingMPIMessages::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::DropIncomingMPIMessages::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::DropIncomingMPIMessages::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::DropIncomingMPIMessages::_log(
    "exahype::mappings::DropIncomingMPIMessages");

#ifdef Parallel
void exahype::mappings::DropIncomingMPIMessages::mergeWithNeighbour(
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
        int receivedMetadataIndex = MetadataHeap::getInstance().createData(0,0);
        assertion(MetadataHeap::getInstance().getData(receivedMetadataIndex).empty());
        MetadataHeap::getInstance().receiveData(
            receivedMetadataIndex,
            fromRank, fineGridX, level,
            peano::heap::MessageType::NeighbourCommunication);
        exahype::MetadataHeap::HeapEntries& receivedMetadata = MetadataHeap::getInstance().getData(receivedMetadataIndex);
        assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());

        if(vertex.hasToMergeWithNeighbourData(src,dest)) {
          vertex.setFaceDataExchangeCountersOfDestination(src,dest,TWO_POWER_D);
          vertex.setMergePerformed(src,dest,true);

          dropNeighbourData(
              fromRank,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              src,dest,
              fineGridX,level,
              receivedMetadata);
        } else {
          dropNeighbourData(
              fromRank,
              vertex.getCellDescriptionsIndex()[srcScalar],
              vertex.getCellDescriptionsIndex()[destScalar],
              src,dest,
              fineGridX,level,
              receivedMetadata);
        }
        // Clean up
        MetadataHeap::getInstance().deleteData(receivedMetadataIndex);
      }
    enddforx
  enddforx
}


void exahype::mappings::DropIncomingMPIMessages::dropNeighbourData(
    const int fromRank,
    const int srcCellDescriptionIndex,
    const int destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS,int>& src,
    const tarch::la::Vector<DIMENSIONS,int>& dest,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int level,
    const exahype::MetadataHeap::HeapEntries& receivedMetadata) {
  logDebug(
      "mergeWithNeighbour(...)", "[drop] " <<
      fromRank << " for vertex x=" << x << ", level=" << level <<
      ", src type=" << multiscalelinkedcell::indexToString(destCellDescriptionIndex) <<
      ", src=" << src << ", dest=" << dest);

  assertion(receivedMetadata.size()==solvers::RegisteredSolvers.size());

  for(auto solver : solvers::RegisteredSolvers) {
    solver->dropNeighbourData(fromRank,src,dest,x,level);
  }
}

void exahype::mappings::DropIncomingMPIMessages::mergeWithMaster(
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



}

//
// All methods below are nop.
//
// ====================================

void exahype::mappings::DropIncomingMPIMessages::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::DropIncomingMPIMessages::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::DropIncomingMPIMessages::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::DropIncomingMPIMessages::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::DropIncomingMPIMessages::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
        const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

bool exahype::mappings::DropIncomingMPIMessages::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return false;
}

void exahype::mappings::DropIncomingMPIMessages::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::DropIncomingMPIMessages::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::DropIncomingMPIMessages::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::DropIncomingMPIMessages::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}
#endif

exahype::mappings::DropIncomingMPIMessages::DropIncomingMPIMessages() {}

exahype::mappings::DropIncomingMPIMessages::~DropIncomingMPIMessages() {}

#if defined(SharedMemoryParallelisation)
exahype::mappings::DropIncomingMPIMessages::DropIncomingMPIMessages(
    const DropIncomingMPIMessages& masterThread) {}

void exahype::mappings::DropIncomingMPIMessages::mergeWithWorkerThread(
    const DropIncomingMPIMessages& workerThread) {}
#endif

void exahype::mappings::DropIncomingMPIMessages::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::DropIncomingMPIMessages::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}


void exahype::mappings::DropIncomingMPIMessages::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::DropIncomingMPIMessages::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::DropIncomingMPIMessages::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::DropIncomingMPIMessages::beginIteration(
    exahype::State& solverState) {
  logTraceInWith1Argument("beginIteration(State)", solverState);

  #ifdef Parallel
  exahype::solvers::ADERDGSolver::Heap::getInstance().finishedToSendSynchronousData();
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().finishedToSendSynchronousData();
  DataHeap::getInstance().finishedToSendSynchronousData();
  MetadataHeap::getInstance().finishedToSendSynchronousData();
  #endif

  logTraceOutWith1Argument("beginIteration(State)", solverState);
}

void exahype::mappings::DropIncomingMPIMessages::endIteration(
    exahype::State& solverState) {}

void exahype::mappings::DropIncomingMPIMessages::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}

void exahype::mappings::DropIncomingMPIMessages::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}
