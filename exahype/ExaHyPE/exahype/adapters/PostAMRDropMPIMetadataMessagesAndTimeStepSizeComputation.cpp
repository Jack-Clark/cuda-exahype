#include "exahype/adapters/PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.h"


peano::CommunicationSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::communicationSpecification() {
  return peano::CommunicationSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::communicationSpecification()
   & exahype::mappings::PreProcessing::communicationSpecification()
   & exahype::mappings::TimeStepSizeComputation::communicationSpecification()
   & exahype::mappings::Sending::communicationSpecification()
   & exahype::mappings::PostProcessing::communicationSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::touchVertexLastTimeSpecification()
   & exahype::mappings::PreProcessing::touchVertexLastTimeSpecification()
   & exahype::mappings::TimeStepSizeComputation::touchVertexLastTimeSpecification()
   & exahype::mappings::Sending::touchVertexLastTimeSpecification()
   & exahype::mappings::PostProcessing::touchVertexLastTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::touchVertexFirstTimeSpecification() { 
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::touchVertexFirstTimeSpecification()
   & exahype::mappings::PreProcessing::touchVertexFirstTimeSpecification()
   & exahype::mappings::TimeStepSizeComputation::touchVertexFirstTimeSpecification()
   & exahype::mappings::Sending::touchVertexFirstTimeSpecification()
   & exahype::mappings::PostProcessing::touchVertexFirstTimeSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::enterCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::enterCellSpecification()
   & exahype::mappings::PreProcessing::enterCellSpecification()
   & exahype::mappings::TimeStepSizeComputation::enterCellSpecification()
   & exahype::mappings::Sending::enterCellSpecification()
   & exahype::mappings::PostProcessing::enterCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::leaveCellSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::leaveCellSpecification()
   & exahype::mappings::PreProcessing::leaveCellSpecification()
   & exahype::mappings::TimeStepSizeComputation::leaveCellSpecification()
   & exahype::mappings::Sending::leaveCellSpecification()
   & exahype::mappings::PostProcessing::leaveCellSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::ascendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::ascendSpecification()
   & exahype::mappings::PreProcessing::ascendSpecification()
   & exahype::mappings::TimeStepSizeComputation::ascendSpecification()
   & exahype::mappings::Sending::ascendSpecification()
   & exahype::mappings::PostProcessing::ascendSpecification()

  ;
}


peano::MappingSpecification   exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::descendSpecification() {
  return peano::MappingSpecification::getMinimalSpecification()
   & exahype::mappings::DropIncomingMPIMetadataMessages::descendSpecification()
   & exahype::mappings::PreProcessing::descendSpecification()
   & exahype::mappings::TimeStepSizeComputation::descendSpecification()
   & exahype::mappings::Sending::descendSpecification()
   & exahype::mappings::PostProcessing::descendSpecification()

  ;
}


exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() {
}


exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::~PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() {
}


#if defined(SharedMemoryParallelisation)
exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation(const PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation&  masterThread):
  _map2DropIncomingMPIMetadataMessages(masterThread._map2DropIncomingMPIMetadataMessages) , 
  _map2PreProcessing(masterThread._map2PreProcessing) , 
  _map2TimeStepSizeComputation(masterThread._map2TimeStepSizeComputation) , 
  _map2Sending(masterThread._map2Sending) , 
  _map2PostProcessing(masterThread._map2PostProcessing) 

{
}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithWorkerThread(const PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation& workerThread) {
  _map2DropIncomingMPIMetadataMessages.mergeWithWorkerThread(workerThread._map2DropIncomingMPIMetadataMessages);
  _map2PreProcessing.mergeWithWorkerThread(workerThread._map2PreProcessing);
  _map2TimeStepSizeComputation.mergeWithWorkerThread(workerThread._map2TimeStepSizeComputation);
  _map2Sending.mergeWithWorkerThread(workerThread._map2Sending);
  _map2PostProcessing.mergeWithWorkerThread(workerThread._map2PostProcessing);

}
#endif


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::createHangingVertex(
      exahype::Vertex&     fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                fineGridH,
      exahype::Vertex * const   coarseGridVertices,
      const peano::grid::VertexEnumerator&      coarseGridVerticesEnumerator,
      exahype::Cell&       coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                   fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );


}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.destroyHangingVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createInnerVertex(fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.createBoundaryVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.destroyVertex( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2DropIncomingMPIMetadataMessages.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PreProcessing.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2TimeStepSizeComputation.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Sending.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.createCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2DropIncomingMPIMetadataMessages.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PreProcessing.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2TimeStepSizeComputation.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Sending.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.destroyCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


#ifdef Parallel
void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithNeighbour(
  exahype::Vertex&  vertex,
  const exahype::Vertex&  neighbour,
  int                                           fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridX,
  const tarch::la::Vector<DIMENSIONS,double>&   fineGridH,
  int                                           level
) {
   _map2DropIncomingMPIMetadataMessages.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2PreProcessing.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2TimeStepSizeComputation.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2Sending.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );
   _map2PostProcessing.mergeWithNeighbour( vertex, neighbour, fromRank, fineGridX, fineGridH, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::prepareSendToNeighbour(
  exahype::Vertex&  vertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2DropIncomingMPIMetadataMessages.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2PreProcessing.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2TimeStepSizeComputation.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2Sending.prepareSendToNeighbour( vertex, toRank, x, h, level );
   _map2PostProcessing.prepareSendToNeighbour( vertex, toRank, x, h, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::prepareCopyToRemoteNode(
  exahype::Vertex&  localVertex,
  int                                           toRank,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2DropIncomingMPIMetadataMessages.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2PreProcessing.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2TimeStepSizeComputation.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2Sending.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );
   _map2PostProcessing.prepareCopyToRemoteNode( localVertex, toRank, x, h, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::prepareCopyToRemoteNode(
  exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
) {
   _map2DropIncomingMPIMetadataMessages.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2PreProcessing.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2TimeStepSizeComputation.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2Sending.prepareCopyToRemoteNode( localCell, toRank, x, h, level );
   _map2PostProcessing.prepareCopyToRemoteNode( localCell, toRank, x, h, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Vertex&  localVertex,
  const exahype::Vertex&  masterOrWorkerVertex,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2DropIncomingMPIMetadataMessages.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2PreProcessing.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2TimeStepSizeComputation.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2Sending.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );
   _map2PostProcessing.mergeWithRemoteDataDueToForkOrJoin( localVertex, masterOrWorkerVertex, fromRank, x, h, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithRemoteDataDueToForkOrJoin(
  exahype::Cell&  localCell,
  const exahype::Cell&  masterOrWorkerCell,
  int                                       fromRank,
  const tarch::la::Vector<DIMENSIONS,double>&  x,
  const tarch::la::Vector<DIMENSIONS,double>&  h,
  int                                       level
) {
   _map2DropIncomingMPIMetadataMessages.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2PreProcessing.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2TimeStepSizeComputation.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2Sending.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );
   _map2PostProcessing.mergeWithRemoteDataDueToForkOrJoin( localCell, masterOrWorkerCell, fromRank, x, h, level );

}


bool exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::prepareSendToWorker(
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
   result |= _map2DropIncomingMPIMetadataMessages.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2PreProcessing.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2TimeStepSizeComputation.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2Sending.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );
   result |= _map2PostProcessing.prepareSendToWorker( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker );

  return result;
}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::prepareSendToMaster(
  exahype::Cell&                       localCell,
  exahype::Vertex *                    vertices,
  const peano::grid::VertexEnumerator&       verticesEnumerator, 
  const exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
  const exahype::Cell&                 coarseGridCell,
  const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
) {
   _map2DropIncomingMPIMetadataMessages.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2PreProcessing.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2TimeStepSizeComputation.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2Sending.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
   _map2PostProcessing.prepareSendToMaster( localCell, vertices, verticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithMaster(
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
   _map2DropIncomingMPIMetadataMessages.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2PreProcessing.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2TimeStepSizeComputation.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2Sending.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );
   _map2PostProcessing.mergeWithMaster( workerGridCell, workerGridVertices, workerEnumerator, fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell, worker, workerState, masterState );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::receiveDataFromMaster(
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
   _map2DropIncomingMPIMetadataMessages.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2PreProcessing.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2TimeStepSizeComputation.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2Sending.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );
   _map2PostProcessing.receiveDataFromMaster( receivedCell, receivedVertices, receivedVerticesEnumerator, receivedCoarseGridVertices, receivedCoarseGridVerticesEnumerator, receivedCoarseGridCell, workersCoarseGridVertices, workersCoarseGridVerticesEnumerator, workersCoarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithWorker(
  exahype::Cell&           localCell, 
  const exahype::Cell&     receivedMasterCell,
  const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
  const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
  int                                          level
) {
   _map2DropIncomingMPIMetadataMessages.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2PreProcessing.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2TimeStepSizeComputation.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2Sending.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );
   _map2PostProcessing.mergeWithWorker( localCell, receivedMasterCell, cellCentre, cellSize, level );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::mergeWithWorker(
  exahype::Vertex&        localVertex,
  const exahype::Vertex&  receivedMasterVertex,
  const tarch::la::Vector<DIMENSIONS,double>&   x,
  const tarch::la::Vector<DIMENSIONS,double>&   h,
  int                                           level
) {
   _map2DropIncomingMPIMetadataMessages.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2PreProcessing.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2TimeStepSizeComputation.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2Sending.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );
   _map2PostProcessing.mergeWithWorker( localVertex, receivedMasterVertex, x, h, level );

}
#endif


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.touchVertexFirstTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
) {
  _map2DropIncomingMPIMetadataMessages.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PreProcessing.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2TimeStepSizeComputation.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2Sending.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );
  _map2PostProcessing.touchVertexLastTime( fineGridVertex, fineGridX, fineGridH, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfVertex );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
) {
  _map2DropIncomingMPIMetadataMessages.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PreProcessing.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2TimeStepSizeComputation.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Sending.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.enterCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::leaveCell(
      exahype::Cell&           fineGridCell,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfCell
) {
  _map2DropIncomingMPIMetadataMessages.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PreProcessing.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2TimeStepSizeComputation.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2Sending.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );
  _map2PostProcessing.leaveCell( fineGridCell, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell, fineGridPositionOfCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::beginIteration(
  exahype::State&  solverState
) {
  _map2DropIncomingMPIMetadataMessages.beginIteration( solverState );
  _map2PreProcessing.beginIteration( solverState );
  _map2TimeStepSizeComputation.beginIteration( solverState );
  _map2Sending.beginIteration( solverState );
  _map2PostProcessing.beginIteration( solverState );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::endIteration(
  exahype::State&  solverState
) {
  _map2DropIncomingMPIMetadataMessages.endIteration( solverState );
  _map2PreProcessing.endIteration( solverState );
  _map2TimeStepSizeComputation.endIteration( solverState );
  _map2Sending.endIteration( solverState );
  _map2PostProcessing.endIteration( solverState );

}




void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::descend(
  exahype::Cell * const          fineGridCells,
  exahype::Vertex * const        fineGridVertices,
  const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
  exahype::Vertex * const        coarseGridVertices,
  const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
  exahype::Cell&                 coarseGridCell
) {
  _map2DropIncomingMPIMetadataMessages.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PreProcessing.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2TimeStepSizeComputation.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Sending.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PostProcessing.descend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}


void exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation::ascend(
  exahype::Cell * const    fineGridCells,
  exahype::Vertex * const  fineGridVertices,
  const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
  exahype::Vertex * const  coarseGridVertices,
  const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
  exahype::Cell&           coarseGridCell
) {
  _map2DropIncomingMPIMetadataMessages.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PreProcessing.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2TimeStepSizeComputation.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2Sending.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );
  _map2PostProcessing.ascend( fineGridCells, fineGridVertices, fineGridVerticesEnumerator, coarseGridVertices, coarseGridVerticesEnumerator, coarseGridCell );

}
