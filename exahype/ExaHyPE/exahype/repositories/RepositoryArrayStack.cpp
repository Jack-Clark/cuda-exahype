#include "exahype/repositories/RepositoryArrayStack.h"

#include "tarch/Assertions.h"
#include "tarch/timing/Watch.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"

#ifdef Parallel
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#endif

#include "peano/datatraversal/autotuning/Oracle.h"

#include "tarch/compiler/CompilerSpecificSettings.h"

#if !defined(CompilerICC)
#include "peano/grid/Grid.cpph"
#endif


tarch::logging::Log exahype::repositories::RepositoryArrayStack::_log( "exahype::repositories::RepositoryArrayStack" );


exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack, maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusMergingAndSpreadingMPI(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusMergingMPI(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithReinitialisation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionRecomputationAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithNeighbourDataMerging(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlot(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlot2d(_vertexStack,_cellStack,_geometry,_solverState,domainSize,domainOffset,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(...)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(...)" );
}



exahype::repositories::RepositoryArrayStack::RepositoryArrayStack(
  peano::geometry::Geometry&                   geometry,
  int                                          maximumSizeOfCellInOutStack,
  int                                          maximumSizeOfVertexInOutStack,
  int                                          maximumSizeOfVertexTemporaryStack
):
  _geometry(geometry),
  _cellStack(maximumSizeOfCellInOutStack),
  _vertexStack(maximumSizeOfVertexInOutStack,maximumSizeOfVertexTemporaryStack),
  _solverState(),
  _gridWithMeshRefinement(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAugmentedAMRGrid(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithInitialConditionAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithGridErasing(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPlotAndADERDGTimeStep(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionRerun(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusSpreading(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusMergingAndSpreadingMPI(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithLimiterStatusMergingMPI(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithReinitialisation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionRecomputationAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithNeighbourDataMerging(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithSolutionUpdate(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithTimeStepSizeComputation(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPrediction(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlot(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),
  _gridWithPredictionAndPlot2d(_vertexStack,_cellStack,_geometry,_solverState,_regularGridContainer,_traversalOrderOnTopLevel),

  _repositoryState() {
  logTraceIn( "RepositoryArrayStack(Geometry&)" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );

  peano::datatraversal::autotuning::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #ifdef Parallel
  peano::parallel::loadbalancing::Oracle::getInstance().setNumberOfOracles(exahype::records::RepositoryState::NumberOfAdapters);
  #endif
  
  logTraceOut( "RepositoryArrayStack(Geometry&)" );
}
    
   
exahype::repositories::RepositoryArrayStack::~RepositoryArrayStack() {
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );
}


void exahype::repositories::RepositoryArrayStack::restart(
  const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
  const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
  int                                          domainLevel,
  const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
) {
  logTraceInWith4Arguments( "restart(...)", domainSize, domainOffset, domainLevel, positionOfCentralElementWithRespectToCoarserRemoteLevel );
  #ifdef Parallel
  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());
  #endif
  
  logInfo( "restart(...)", "start node for subdomain " << domainOffset << "x" << domainSize << " on level " << domainLevel << " with master " << tarch::parallel::NodePool::getInstance().getMasterRank() );
  
  assertion( _repositoryState.getAction() == exahype::records::RepositoryState::Terminate );

  _vertexStack.clear();
  _cellStack.clear();

  _gridWithMeshRefinement.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAugmentedAMRGrid.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithInitialConditionAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndFusedTimeSteppingInitialisation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithGridErasing.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithADERDGTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPlotAndADERDGTimeStep.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionRerun.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithLimiterStatusSpreading.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithLimiterStatusMergingAndSpreadingMPI.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithLimiterStatusMergingMPI.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithReinitialisation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionRecomputationAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithNeighbourDataMerging.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithSolutionUpdate.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithTimeStepSizeComputation.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPrediction.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndPlot.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);
  _gridWithPredictionAndPlot2d.restart(domainSize,domainOffset,domainLevel,positionOfCentralElementWithRespectToCoarserRemoteLevel);

 
   _solverState.restart();
 
  logTraceOut( "restart(...)" );
}


void exahype::repositories::RepositoryArrayStack::terminate() {
  logTraceIn( "terminate()" );
  
  _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  peano::parallel::SendReceiveBufferPool::getInstance().terminate();
  #endif
  
  _gridWithMeshRefinement.terminate();
  _gridWithPlotAugmentedAMRGrid.terminate();
  _gridWithInitialConditionAndTimeStepSizeComputation.terminate();
  _gridWithPredictionAndFusedTimeSteppingInitialisation.terminate();
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot.terminate();
  _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d.terminate();
  _gridWithGridErasing.terminate();
  _gridWithADERDGTimeStep.terminate();
  _gridWithPlotAndADERDGTimeStep.terminate();
  _gridWithPredictionRerun.terminate();
  _gridWithLimiterStatusSpreading.terminate();
  _gridWithLimiterStatusMergingAndSpreadingMPI.terminate();
  _gridWithLimiterStatusMergingMPI.terminate();
  _gridWithReinitialisation.terminate();
  _gridWithSolutionRecomputationAndTimeStepSizeComputation.terminate();
  _gridWithNeighbourDataMerging.terminate();
  _gridWithSolutionUpdate.terminate();
  _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.terminate();
  _gridWithTimeStepSizeComputation.terminate();
  _gridWithPrediction.terminate();
  _gridWithPredictionAndPlot.terminate();
  _gridWithPredictionAndPlot2d.terminate();

 
  logTraceOut( "terminate()" );
}


exahype::State& exahype::repositories::RepositoryArrayStack::getState() {
  return _solverState;
}


const exahype::State& exahype::repositories::RepositoryArrayStack::getState() const {
  return _solverState;
}

   
void exahype::repositories::RepositoryArrayStack::iterate(int numberOfIterations, bool exchangeBoundaryVertices) {
  tarch::timing::Watch watch( "exahype::repositories::RepositoryArrayStack", "iterate(bool)", false);
  
  #ifdef Parallel
  if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
    _repositoryState.setNumberOfIterations(numberOfIterations);
    _repositoryState.setExchangeBoundaryVertices(exchangeBoundaryVertices);
    tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
      _repositoryState,
      peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
    );
  }
  else {
    assertionEquals( numberOfIterations, 1 );
    numberOfIterations = _repositoryState.getNumberOfIterations();
  }

  peano::parallel::SendReceiveBufferPool::getInstance().exchangeBoundaryVertices(_repositoryState.getExchangeBoundaryVertices());

  if ( numberOfIterations > 1 && _solverState.isInvolvedInJoinOrFork() ) {
    logWarning( "iterate()", "iterate invoked for multiple traversals though load balancing still does redistribute data" );
  }
  bool switchedLoadBalancingTemporarilyOff = false;
  if ( numberOfIterations > 1 && peano::parallel::loadbalancing::Oracle::getInstance().isLoadBalancingActivated() ) {
    switchedLoadBalancingTemporarilyOff = true;
    peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(false);
  }

  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());

  peano::parallel::loadbalancing::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  
  _solverState.currentlyRunsMultipleIterations(_repositoryState.getNumberOfIterations()>1);
  #else
  peano::datatraversal::autotuning::Oracle::getInstance().switchToOracle(_repositoryState.getAction());
  #endif
  
  for (int i=0; i<numberOfIterations; i++) {
    switch ( _repositoryState.getAction()) {
      case exahype::records::RepositoryState::UseAdapterMeshRefinement: watch.startTimer(); _gridWithMeshRefinement.iterate(); watch.stopTimer(); _measureMeshRefinementCPUTime.setValue( watch.getCPUTime() ); _measureMeshRefinementCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid: watch.startTimer(); _gridWithPlotAugmentedAMRGrid.iterate(); watch.stopTimer(); _measurePlotAugmentedAMRGridCPUTime.setValue( watch.getCPUTime() ); _measurePlotAugmentedAMRGridCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation: watch.startTimer(); _gridWithInitialConditionAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureInitialConditionAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureInitialConditionAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisation: watch.startTimer(); _gridWithPredictionAndFusedTimeSteppingInitialisation.iterate(); watch.stopTimer(); _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot: watch.startTimer(); _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot.iterate(); watch.stopTimer(); _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d: watch.startTimer(); _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d.iterate(); watch.stopTimer(); _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterGridErasing: watch.startTimer(); _gridWithGridErasing.iterate(); watch.stopTimer(); _measureGridErasingCPUTime.setValue( watch.getCPUTime() ); _measureGridErasingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterADERDGTimeStep: watch.startTimer(); _gridWithADERDGTimeStep.iterate(); watch.stopTimer(); _measureADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measureADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep: watch.startTimer(); _gridWithPlotAndADERDGTimeStep.iterate(); watch.stopTimer(); _measurePlotAndADERDGTimeStepCPUTime.setValue( watch.getCPUTime() ); _measurePlotAndADERDGTimeStepCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionRerun: watch.startTimer(); _gridWithPredictionRerun.iterate(); watch.stopTimer(); _measurePredictionRerunCPUTime.setValue( watch.getCPUTime() ); _measurePredictionRerunCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading: watch.startTimer(); _gridWithLimiterStatusSpreading.iterate(); watch.stopTimer(); _measureLimiterStatusSpreadingCPUTime.setValue( watch.getCPUTime() ); _measureLimiterStatusSpreadingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterLimiterStatusMergingAndSpreadingMPI: watch.startTimer(); _gridWithLimiterStatusMergingAndSpreadingMPI.iterate(); watch.stopTimer(); _measureLimiterStatusMergingAndSpreadingMPICPUTime.setValue( watch.getCPUTime() ); _measureLimiterStatusMergingAndSpreadingMPICalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterLimiterStatusMergingMPI: watch.startTimer(); _gridWithLimiterStatusMergingMPI.iterate(); watch.stopTimer(); _measureLimiterStatusMergingMPICPUTime.setValue( watch.getCPUTime() ); _measureLimiterStatusMergingMPICalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterReinitialisation: watch.startTimer(); _gridWithReinitialisation.iterate(); watch.stopTimer(); _measureReinitialisationCPUTime.setValue( watch.getCPUTime() ); _measureReinitialisationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation: watch.startTimer(); _gridWithSolutionRecomputationAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterNeighbourDataMerging: watch.startTimer(); _gridWithNeighbourDataMerging.iterate(); watch.stopTimer(); _measureNeighbourDataMergingCPUTime.setValue( watch.getCPUTime() ); _measureNeighbourDataMergingCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterSolutionUpdate: watch.startTimer(); _gridWithSolutionUpdate.iterate(); watch.stopTimer(); _measureSolutionUpdateCPUTime.setValue( watch.getCPUTime() ); _measureSolutionUpdateCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation: watch.startTimer(); _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.iterate(); watch.stopTimer(); _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation: watch.startTimer(); _gridWithTimeStepSizeComputation.iterate(); watch.stopTimer(); _measureTimeStepSizeComputationCPUTime.setValue( watch.getCPUTime() ); _measureTimeStepSizeComputationCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPrediction: watch.startTimer(); _gridWithPrediction.iterate(); watch.stopTimer(); _measurePredictionCPUTime.setValue( watch.getCPUTime() ); _measurePredictionCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndPlot: watch.startTimer(); _gridWithPredictionAndPlot.iterate(); watch.stopTimer(); _measurePredictionAndPlotCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndPlotCalendarTime.setValue( watch.getCalendarTime() ); break;
      case exahype::records::RepositoryState::UseAdapterPredictionAndPlot2d: watch.startTimer(); _gridWithPredictionAndPlot2d.iterate(); watch.stopTimer(); _measurePredictionAndPlot2dCPUTime.setValue( watch.getCPUTime() ); _measurePredictionAndPlot2dCalendarTime.setValue( watch.getCalendarTime() ); break;

      case exahype::records::RepositoryState::Terminate:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::NumberOfAdapters:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::RunOnAllNodes:
        assertionMsg( false, "this branch/state should never be reached" ); 
        break;
      case exahype::records::RepositoryState::ReadCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
      case exahype::records::RepositoryState::WriteCheckpoint:
        assertionMsg( false, "not implemented yet" );
        break;
    }
    #ifdef Parallel
    if ( switchedLoadBalancingTemporarilyOff && i==numberOfIterations-1) {
      peano::parallel::loadbalancing::Oracle::getInstance().activateLoadBalancing(true);
    }
    #endif
  }
    
  #ifdef Parallel
  if (_solverState.isJoiningWithMaster()) {
    _repositoryState.setAction( exahype::records::RepositoryState::Terminate );
  }
  #endif
}

 void exahype::repositories::RepositoryArrayStack::switchToMeshRefinement() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterMeshRefinement); }
 void exahype::repositories::RepositoryArrayStack::switchToPlotAugmentedAMRGrid() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid); }
 void exahype::repositories::RepositoryArrayStack::switchToInitialConditionAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndFusedTimeSteppingInitialisation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisation); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndFusedTimeSteppingInitialisationAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndFusedTimeSteppingInitialisationAndPlot2d() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d); }
 void exahype::repositories::RepositoryArrayStack::switchToGridErasing() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterGridErasing); }
 void exahype::repositories::RepositoryArrayStack::switchToADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterADERDGTimeStep); }
 void exahype::repositories::RepositoryArrayStack::switchToPlotAndADERDGTimeStep() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionRerun() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionRerun); }
 void exahype::repositories::RepositoryArrayStack::switchToLimiterStatusSpreading() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading); }
 void exahype::repositories::RepositoryArrayStack::switchToLimiterStatusMergingAndSpreadingMPI() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterLimiterStatusMergingAndSpreadingMPI); }
 void exahype::repositories::RepositoryArrayStack::switchToLimiterStatusMergingMPI() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterLimiterStatusMergingMPI); }
 void exahype::repositories::RepositoryArrayStack::switchToReinitialisation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterReinitialisation); }
 void exahype::repositories::RepositoryArrayStack::switchToSolutionRecomputationAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToNeighbourDataMerging() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterNeighbourDataMerging); }
 void exahype::repositories::RepositoryArrayStack::switchToSolutionUpdate() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterSolutionUpdate); }
 void exahype::repositories::RepositoryArrayStack::switchToPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToTimeStepSizeComputation() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation); }
 void exahype::repositories::RepositoryArrayStack::switchToPrediction() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPrediction); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndPlot() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndPlot); }
 void exahype::repositories::RepositoryArrayStack::switchToPredictionAndPlot2d() { _repositoryState.setAction(exahype::records::RepositoryState::UseAdapterPredictionAndPlot2d); }



 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterMeshRefinement() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterMeshRefinement; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlotAugmentedAMRGrid() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAugmentedAMRGrid; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterInitialConditionAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterInitialConditionAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndFusedTimeSteppingInitialisation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterGridErasing() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterGridErasing; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterADERDGTimeStep; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPlotAndADERDGTimeStep() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPlotAndADERDGTimeStep; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionRerun() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionRerun; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterLimiterStatusSpreading() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterLimiterStatusSpreading; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterLimiterStatusMergingAndSpreadingMPI() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterLimiterStatusMergingAndSpreadingMPI; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterLimiterStatusMergingMPI() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterLimiterStatusMergingMPI; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterReinitialisation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterReinitialisation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterSolutionRecomputationAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionRecomputationAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterNeighbourDataMerging() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterNeighbourDataMerging; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterSolutionUpdate() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterSolutionUpdate; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterTimeStepSizeComputation() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterTimeStepSizeComputation; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPrediction() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPrediction; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndPlot() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndPlot; }
 bool exahype::repositories::RepositoryArrayStack::isActiveAdapterPredictionAndPlot2d() const { return _repositoryState.getAction() == exahype::records::RepositoryState::UseAdapterPredictionAndPlot2d; }



peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* exahype::repositories::RepositoryArrayStack::createEmptyCheckpoint() {
  return new peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>();
} 


void exahype::repositories::RepositoryArrayStack::writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint) {
  _solverState.writeToCheckpoint( *checkpoint );
  _vertexStack.writeToCheckpoint( *checkpoint );
  _cellStack.writeToCheckpoint( *checkpoint );
} 


void exahype::repositories::RepositoryArrayStack::setMaximumMemoryFootprintForTemporaryRegularGrids(double value) {
  _regularGridContainer.setMaximumMemoryFootprintForTemporaryRegularGrids(value);
}


void exahype::repositories::RepositoryArrayStack::readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint ) {
  assertionMsg( checkpoint->isValid(), "checkpoint has to be valid if you call this operation" );

  _solverState.readFromCheckpoint( *checkpoint );
  _vertexStack.readFromCheckpoint( *checkpoint );
  _cellStack.readFromCheckpoint( *checkpoint );
}


#ifdef Parallel
void exahype::repositories::RepositoryArrayStack::runGlobalStep() {
  assertion(tarch::parallel::Node::getInstance().isGlobalMaster());

  exahype::records::RepositoryState intermediateStateForWorkingNodes;
  intermediateStateForWorkingNodes.setAction( exahype::records::RepositoryState::RunOnAllNodes );
  
  tarch::parallel::NodePool::getInstance().broadcastToWorkingNodes(
    intermediateStateForWorkingNodes,
    peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag()
  );
  tarch::parallel::NodePool::getInstance().activateIdleNodes();
}


exahype::repositories::RepositoryArrayStack::ContinueCommand exahype::repositories::RepositoryArrayStack::continueToIterate() {
  logTraceIn( "continueToIterate()" );

  assertion( !tarch::parallel::Node::getInstance().isGlobalMaster());

  ContinueCommand result;
  if ( _solverState.hasJoinedWithMaster() ) {
    result = Terminate;
  }
  else {
    int masterNode = tarch::parallel::Node::getInstance().getGlobalMasterRank();
    assertion( masterNode != -1 );

    _repositoryState.receive( masterNode, peano::parallel::SendReceiveBufferPool::getInstance().getIterationManagementTag(), true, ReceiveIterationControlMessagesBlocking );

    result = Continue;
    if (_repositoryState.getAction()==exahype::records::RepositoryState::Terminate) {
      result = Terminate;
    } 
    if (_repositoryState.getAction()==exahype::records::RepositoryState::RunOnAllNodes) {
      result = RunGlobalStep;
    } 
  }
   
  logTraceOutWith1Argument( "continueToIterate()", result );
  return result;
}
#endif


void exahype::repositories::RepositoryArrayStack::logIterationStatistics(bool logAllAdapters) const {
  logInfo( "logIterationStatistics()", "|| adapter name \t || iterations \t || total CPU time [t]=s \t || average CPU time [t]=s \t || total user time [t]=s \t || average user time [t]=s  || CPU time properties  || user time properties " );  
   if (logAllAdapters || _measureMeshRefinementCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| MeshRefinement \t |  " << _measureMeshRefinementCPUTime.getNumberOfMeasurements() << " \t |  " << _measureMeshRefinementCPUTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCPUTime.getValue()  << " \t |  " << _measureMeshRefinementCalendarTime.getAccumulatedValue() << " \t |  " << _measureMeshRefinementCalendarTime.getValue() << " \t |  " << _measureMeshRefinementCPUTime.toString() << " \t |  " << _measureMeshRefinementCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAugmentedAMRGrid \t |  " << _measurePlotAugmentedAMRGridCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.getValue()  << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.getValue() << " \t |  " << _measurePlotAugmentedAMRGridCPUTime.toString() << " \t |  " << _measurePlotAugmentedAMRGridCalendarTime.toString() );
   if (logAllAdapters || _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| InitialConditionAndTimeStepSizeComputation \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureInitialConditionAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndFusedTimeSteppingInitialisation \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.getValue()  << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime.getValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.toString() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndFusedTimeSteppingInitialisationAndPlot \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.getValue()  << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime.getValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.toString() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndFusedTimeSteppingInitialisationAndPlot2d \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.getValue()  << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime.getValue() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.toString() << " \t |  " << _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime.toString() );
   if (logAllAdapters || _measureGridErasingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| GridErasing \t |  " << _measureGridErasingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureGridErasingCPUTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCPUTime.getValue()  << " \t |  " << _measureGridErasingCalendarTime.getAccumulatedValue() << " \t |  " << _measureGridErasingCalendarTime.getValue() << " \t |  " << _measureGridErasingCPUTime.toString() << " \t |  " << _measureGridErasingCalendarTime.toString() );
   if (logAllAdapters || _measureADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| ADERDGTimeStep \t |  " << _measureADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measureADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measureADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measureADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measureADERDGTimeStepCPUTime.toString() << " \t |  " << _measureADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measurePlotAndADERDGTimeStepCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PlotAndADERDGTimeStep \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getAccumulatedValue() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.getValue()  << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.getAccumulatedValue() << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.getValue() << " \t |  " << _measurePlotAndADERDGTimeStepCPUTime.toString() << " \t |  " << _measurePlotAndADERDGTimeStepCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionRerunCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionRerun \t |  " << _measurePredictionRerunCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionRerunCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCPUTime.getValue()  << " \t |  " << _measurePredictionRerunCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionRerunCalendarTime.getValue() << " \t |  " << _measurePredictionRerunCPUTime.toString() << " \t |  " << _measurePredictionRerunCalendarTime.toString() );
   if (logAllAdapters || _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| LimiterStatusSpreading \t |  " << _measureLimiterStatusSpreadingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.getValue()  << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.getValue() << " \t |  " << _measureLimiterStatusSpreadingCPUTime.toString() << " \t |  " << _measureLimiterStatusSpreadingCalendarTime.toString() );
   if (logAllAdapters || _measureLimiterStatusMergingAndSpreadingMPICPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| LimiterStatusMergingAndSpreadingMPI \t |  " << _measureLimiterStatusMergingAndSpreadingMPICPUTime.getNumberOfMeasurements() << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICPUTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICPUTime.getValue()  << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICalendarTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICalendarTime.getValue() << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICPUTime.toString() << " \t |  " << _measureLimiterStatusMergingAndSpreadingMPICalendarTime.toString() );
   if (logAllAdapters || _measureLimiterStatusMergingMPICPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| LimiterStatusMergingMPI \t |  " << _measureLimiterStatusMergingMPICPUTime.getNumberOfMeasurements() << " \t |  " << _measureLimiterStatusMergingMPICPUTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusMergingMPICPUTime.getValue()  << " \t |  " << _measureLimiterStatusMergingMPICalendarTime.getAccumulatedValue() << " \t |  " << _measureLimiterStatusMergingMPICalendarTime.getValue() << " \t |  " << _measureLimiterStatusMergingMPICPUTime.toString() << " \t |  " << _measureLimiterStatusMergingMPICalendarTime.toString() );
   if (logAllAdapters || _measureReinitialisationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Reinitialisation \t |  " << _measureReinitialisationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureReinitialisationCPUTime.getAccumulatedValue() << " \t |  " << _measureReinitialisationCPUTime.getValue()  << " \t |  " << _measureReinitialisationCalendarTime.getAccumulatedValue() << " \t |  " << _measureReinitialisationCalendarTime.getValue() << " \t |  " << _measureReinitialisationCPUTime.toString() << " \t |  " << _measureReinitialisationCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionRecomputationAndTimeStepSizeComputation \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureNeighbourDataMergingCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| NeighbourDataMerging \t |  " << _measureNeighbourDataMergingCPUTime.getNumberOfMeasurements() << " \t |  " << _measureNeighbourDataMergingCPUTime.getAccumulatedValue() << " \t |  " << _measureNeighbourDataMergingCPUTime.getValue()  << " \t |  " << _measureNeighbourDataMergingCalendarTime.getAccumulatedValue() << " \t |  " << _measureNeighbourDataMergingCalendarTime.getValue() << " \t |  " << _measureNeighbourDataMergingCPUTime.toString() << " \t |  " << _measureNeighbourDataMergingCalendarTime.toString() );
   if (logAllAdapters || _measureSolutionUpdateCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| SolutionUpdate \t |  " << _measureSolutionUpdateCPUTime.getNumberOfMeasurements() << " \t |  " << _measureSolutionUpdateCPUTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCPUTime.getValue()  << " \t |  " << _measureSolutionUpdateCalendarTime.getAccumulatedValue() << " \t |  " << _measureSolutionUpdateCalendarTime.getValue() << " \t |  " << _measureSolutionUpdateCPUTime.toString() << " \t |  " << _measureSolutionUpdateCalendarTime.toString() );
   if (logAllAdapters || _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measureTimeStepSizeComputationCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| TimeStepSizeComputation \t |  " << _measureTimeStepSizeComputationCPUTime.getNumberOfMeasurements() << " \t |  " << _measureTimeStepSizeComputationCPUTime.getAccumulatedValue() << " \t |  " << _measureTimeStepSizeComputationCPUTime.getValue()  << " \t |  " << _measureTimeStepSizeComputationCalendarTime.getAccumulatedValue() << " \t |  " << _measureTimeStepSizeComputationCalendarTime.getValue() << " \t |  " << _measureTimeStepSizeComputationCPUTime.toString() << " \t |  " << _measureTimeStepSizeComputationCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| Prediction \t |  " << _measurePredictionCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionCPUTime.getValue()  << " \t |  " << _measurePredictionCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionCalendarTime.getValue() << " \t |  " << _measurePredictionCPUTime.toString() << " \t |  " << _measurePredictionCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndPlotCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndPlot \t |  " << _measurePredictionAndPlotCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndPlotCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotCPUTime.getValue()  << " \t |  " << _measurePredictionAndPlotCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlotCalendarTime.getValue() << " \t |  " << _measurePredictionAndPlotCPUTime.toString() << " \t |  " << _measurePredictionAndPlotCalendarTime.toString() );
   if (logAllAdapters || _measurePredictionAndPlot2dCPUTime.getNumberOfMeasurements()>0) logInfo( "logIterationStatistics()", "| PredictionAndPlot2d \t |  " << _measurePredictionAndPlot2dCPUTime.getNumberOfMeasurements() << " \t |  " << _measurePredictionAndPlot2dCPUTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlot2dCPUTime.getValue()  << " \t |  " << _measurePredictionAndPlot2dCalendarTime.getAccumulatedValue() << " \t |  " << _measurePredictionAndPlot2dCalendarTime.getValue() << " \t |  " << _measurePredictionAndPlot2dCPUTime.toString() << " \t |  " << _measurePredictionAndPlot2dCalendarTime.toString() );

}


void exahype::repositories::RepositoryArrayStack::clearIterationStatistics() {
   _measureMeshRefinementCPUTime.erase();
   _measurePlotAugmentedAMRGridCPUTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCPUTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationCPUTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime.erase();
   _measureGridErasingCPUTime.erase();
   _measureADERDGTimeStepCPUTime.erase();
   _measurePlotAndADERDGTimeStepCPUTime.erase();
   _measurePredictionRerunCPUTime.erase();
   _measureLimiterStatusSpreadingCPUTime.erase();
   _measureLimiterStatusMergingAndSpreadingMPICPUTime.erase();
   _measureLimiterStatusMergingMPICPUTime.erase();
   _measureReinitialisationCPUTime.erase();
   _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime.erase();
   _measureNeighbourDataMergingCPUTime.erase();
   _measureSolutionUpdateCPUTime.erase();
   _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime.erase();
   _measureTimeStepSizeComputationCPUTime.erase();
   _measurePredictionCPUTime.erase();
   _measurePredictionAndPlotCPUTime.erase();
   _measurePredictionAndPlot2dCPUTime.erase();

   _measureMeshRefinementCalendarTime.erase();
   _measurePlotAugmentedAMRGridCalendarTime.erase();
   _measureInitialConditionAndTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime.erase();
   _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime.erase();
   _measureGridErasingCalendarTime.erase();
   _measureADERDGTimeStepCalendarTime.erase();
   _measurePlotAndADERDGTimeStepCalendarTime.erase();
   _measurePredictionRerunCalendarTime.erase();
   _measureLimiterStatusSpreadingCalendarTime.erase();
   _measureLimiterStatusMergingAndSpreadingMPICalendarTime.erase();
   _measureLimiterStatusMergingMPICalendarTime.erase();
   _measureReinitialisationCalendarTime.erase();
   _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime.erase();
   _measureNeighbourDataMergingCalendarTime.erase();
   _measureSolutionUpdateCalendarTime.erase();
   _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime.erase();
   _measureTimeStepSizeComputationCalendarTime.erase();
   _measurePredictionCalendarTime.erase();
   _measurePredictionAndPlotCalendarTime.erase();
   _measurePredictionAndPlot2dCalendarTime.erase();

}
