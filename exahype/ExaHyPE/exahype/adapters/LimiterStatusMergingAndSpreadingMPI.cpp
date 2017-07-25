#include "exahype/adapters/LimiterStatusMergingAndSpreadingMPI.h"


peano::CommunicationSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::communicationSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexLastTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::touchVertexFirstTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::enterCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::leaveCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::ascendSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::LimiterStatusMergingAndSpreadingMPI::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::LimiterStatusMergingAndSpreadingMPI::descendSpecification()

  ;
}


exahype::adapters::LimiterStatusMergingAndSpreadingMPI::LimiterStatusMergingAndSpreadingMPI() {
}


exahype::adapters::LimiterStatusMergingAndSpreadingMPI::~LimiterStatusMergingAndSpreadingMPI() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::LimiterStatusMergingAndSpreadingMPI::LimiterStatusMergingAndSpreadingMPI(const LimiterStatusMergingAndSpreadingMPI&  masterThread):
  _map2LimiterStatusMergingAndSpreadingMPI(masterThread._map2LimiterStatusMergingAndSpreadingMPI) 

{
}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithWorkerThread(const LimiterStatusMergingAndSpreadingMPI& workerThread) {
  _map2LimiterStatusMergingAndSpreadingMPI.mergeWithWorkerThread(workerThread._map2LimiterStatusMergingAndSpreadingMPI);

}
#endif


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


#ifdef Parallel
void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.prepareSendToNeighbour( vertex, toRank, x, h, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}


bool exahype::adapters::LimiterStatusMergingAndSpreadingMPI::prepareSendToWorker(
  exahype::Cell&                 fineGridCell,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
  int                                                                  worker
) {
  bool result = false;
   result |= _map2LimiterStatusMergingAndSpreadingMPI.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
   _map2LimiterStatusMergingAndSpreadingMPI.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithMaster(
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
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::receiveDataFromMaster(
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
   _map2LimiterStatusMergingAndSpreadingMPI.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2LimiterStatusMergingAndSpreadingMPI.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2LimiterStatusMergingAndSpreadingMPI.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::beginIteration(
  exahype::State&  solverState
) {
  _map2LimiterStatusMergingAndSpreadingMPI.beginIteration( solverState );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::endIteration(
  exahype::State&  solverState
) {
  _map2LimiterStatusMergingAndSpreadingMPI.endIteration( solverState );

}




void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void exahype::adapters::LimiterStatusMergingAndSpreadingMPI::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  _map2LimiterStatusMergingAndSpreadingMPI.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
