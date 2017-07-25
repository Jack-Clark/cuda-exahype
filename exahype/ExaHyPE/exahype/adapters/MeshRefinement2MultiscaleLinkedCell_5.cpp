#include "exahype/adapters/MeshRefinement2MultiscaleLinkedCell_5.h"

#include <sstream>

#include "peano/utils/Loop.h"
#include "peano/grid/CellFlags.h"


#include "multiscalelinkedcell/HangingVertexBookkeeper.h"
#include "exahype/VertexOperations.h"


peano::CommunicationSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification();
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::leaveCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::ascendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


peano::MappingSpecification   exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::descendSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::Nop,peano::MappingSpecification::AvoidFineGridRaces,false);
}


exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::MeshRefinement2MultiscaleLinkedCell_5() {
}


exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::~MeshRefinement2MultiscaleLinkedCell_5() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::MeshRefinement2MultiscaleLinkedCell_5(const MeshRefinement2MultiscaleLinkedCell_5&  masterThread) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithWorkerThread(const MeshRefinement2MultiscaleLinkedCell_5& workerThread) {
}
#endif


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::createHangingVertex(
  exahype::Vertex&     fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
  exahype::Vertex * const   coarseGridVertices,
  const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
  exahype::Cell&       coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  const int level = coarseGridVerticesEnumerator.getLevel()+1;
  
  fineGridVertex.getCellDescriptionsIndex() = 
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createHangingVertex(
      fineGridX,level,
      fineGridPositionOfVertex,
      VertexOperations::readCellDescriptionsIndex(coarseGridVerticesEnumerator,coarseGridVertices)
    );
}



void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::destroyHangingVertex(
  const exahype::Vertex&   fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::createInnerVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  fineGridVertex.getCellDescriptionsIndex() = multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForNewVertex();
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::createBoundaryVertex(
  exahype::Vertex&               fineGridVertex,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  fineGridVertex.getCellDescriptionsIndex() = multiscalelinkedcell::HangingVertexBookkeeper::getInstance().createVertexLinkMapForBoundaryVertex();
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::createCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::destroyCell(
  const exahype::Cell&           fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().destroyCell(fineGridCell.getCellDescriptionsIndex());
}


#ifdef Parallel
void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
  VertexOperations::writeCellDescriptionsIndex(
    vertex,
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
      vertex.getAdjacentRanks(),
      VertexOperations::readCellDescriptionsIndex(vertex)
    )
  );
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
}


bool exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  return false;
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithMaster(
  const exahype::Cell&           workerGridCell,
  exahype::Vertex * const        workerGridVertices,
  const peano::grid::VertexEnumerator& workerEnumerator,
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker,
  const exahype::State&          workerState,
  exahype::State&                masterState
) {
  dfor2(k)
    VertexOperations::writeCellDescriptionsIndex(
      fineGridVertices[ fineGridVerticesEnumerator(k) ],
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
        fineGridVertices[ fineGridVerticesEnumerator(k) ].getAdjacentRanks(),
        VertexOperations::readCellDescriptionsIndex(fineGridVertices[ fineGridVerticesEnumerator(k) ])
      )
    );
  enddforx
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::receiveDataFromMaster(
      exahype::Cell&                        receivedCell, 
      exahype::Vertex *                     receivedVertices,
      const peano::grid::VertexEnumerator&        receivedVerticesEnumerator,
      exahype::Vertex * const               receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator&        receivedCoarseGridVerticesEnumerator,
      exahype::Cell&                        receivedCoarseGridCell,
      exahype::Vertex * const               workersCoarseGridVertices,
      const peano::grid::VertexEnumerator&        workersCoarseGridVerticesEnumerator,
      exahype::Cell&                        workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&    fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
  VertexOperations::writeCellDescriptionsIndex(
    localVertex,
    multiscalelinkedcell::HangingVertexBookkeeper::getInstance().updateCellIndicesInMergeWithNeighbour(
      localVertex.getAdjacentRanks(),
      VertexOperations::readCellDescriptionsIndex(localVertex)
    )
  );
}
#endif


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::enterCell(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  dfor2(k)
    if (fineGridVertices[fineGridVerticesEnumerator(k)].isHangingNode()) {
      multiscalelinkedcell::HangingVertexBookkeeper::getInstance().getAdjacencyEntriesOfVertex( 
        fineGridVerticesEnumerator.getVertexPosition(k),
        fineGridVerticesEnumerator.getLevel()
      )(TWO_POWER_D-kScalar-1) = fineGridCell.getCellDescriptionsIndex();
    }
    else {
      fineGridVertices[fineGridVerticesEnumerator(k)].getCellDescriptionsIndex()(TWO_POWER_D-kScalar-1) = fineGridCell.getCellDescriptionsIndex();
    }
  enddforx
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::beginIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().beginIteration();
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::endIteration(
  exahype::State&  solverState
) {
  multiscalelinkedcell::HangingVertexBookkeeper::getInstance().endIteration();
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
}


void exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
}
