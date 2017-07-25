// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 
#define _EXAHYPE_REPOSITORIES_REPOSITORY_ARRAY_STD_H_ 


#include "exahype/repositories/Repository.h"
#include "exahype/records/RepositoryState.h"

#include "exahype/State.h"
#include "exahype/Vertex.h"
#include "exahype/Cell.h"

#include "peano/grid/Grid.h"
#include "peano/stacks/CellSTDStack.h"
#include "peano/stacks/VertexSTDStack.h"


 #include "exahype/adapters/MeshRefinement.h" 
 #include "exahype/adapters/PlotAugmentedAMRGrid.h" 
 #include "exahype/adapters/InitialConditionAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisation.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisationAndPlot.h" 
 #include "exahype/adapters/PredictionAndFusedTimeSteppingInitialisationAndPlot2d.h" 
 #include "exahype/adapters/GridErasing.h" 
 #include "exahype/adapters/ADERDGTimeStep.h" 
 #include "exahype/adapters/PlotAndADERDGTimeStep.h" 
 #include "exahype/adapters/PredictionRerun.h" 
 #include "exahype/adapters/LimiterStatusSpreading.h" 
 #include "exahype/adapters/LimiterStatusMergingAndSpreadingMPI.h" 
 #include "exahype/adapters/LimiterStatusMergingMPI.h" 
 #include "exahype/adapters/Reinitialisation.h" 
 #include "exahype/adapters/SolutionRecomputationAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/NeighbourDataMerging.h" 
 #include "exahype/adapters/SolutionUpdate.h" 
 #include "exahype/adapters/PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation.h" 
 #include "exahype/adapters/TimeStepSizeComputation.h" 
 #include "exahype/adapters/Prediction.h" 
 #include "exahype/adapters/PredictionAndPlot.h" 
 #include "exahype/adapters/PredictionAndPlot2d.h" 



namespace exahype {
      namespace repositories {
        class RepositorySTDStack;  
      }
}


class exahype::repositories::RepositorySTDStack: public exahype::repositories::Repository {
  private:
    static tarch::logging::Log _log;
  
    peano::geometry::Geometry& _geometry;
    
    typedef peano::stacks::CellSTDStack<exahype::Cell>       CellStack;
    typedef peano::stacks::VertexSTDStack<exahype::Vertex>   VertexStack;

    CellStack    _cellStack;
    VertexStack  _vertexStack;
    exahype::State          _solverState;
    peano::grid::RegularGridContainer<exahype::Vertex,exahype::Cell>  _regularGridContainer;
    peano::grid::TraversalOrderOnTopLevel                                         _traversalOrderOnTopLevel;

    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::MeshRefinement> _gridWithMeshRefinement;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAugmentedAMRGrid> _gridWithPlotAugmentedAMRGrid;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::InitialConditionAndTimeStepSizeComputation> _gridWithInitialConditionAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndFusedTimeSteppingInitialisation> _gridWithPredictionAndFusedTimeSteppingInitialisation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndFusedTimeSteppingInitialisationAndPlot> _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndFusedTimeSteppingInitialisationAndPlot2d> _gridWithPredictionAndFusedTimeSteppingInitialisationAndPlot2d;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::GridErasing> _gridWithGridErasing;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::ADERDGTimeStep> _gridWithADERDGTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PlotAndADERDGTimeStep> _gridWithPlotAndADERDGTimeStep;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionRerun> _gridWithPredictionRerun;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusSpreading> _gridWithLimiterStatusSpreading;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusMergingAndSpreadingMPI> _gridWithLimiterStatusMergingAndSpreadingMPI;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::LimiterStatusMergingMPI> _gridWithLimiterStatusMergingMPI;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Reinitialisation> _gridWithReinitialisation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionRecomputationAndTimeStepSizeComputation> _gridWithSolutionRecomputationAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::NeighbourDataMerging> _gridWithNeighbourDataMerging;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::SolutionUpdate> _gridWithSolutionUpdate;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation> _gridWithPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::TimeStepSizeComputation> _gridWithTimeStepSizeComputation;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::Prediction> _gridWithPrediction;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndPlot> _gridWithPredictionAndPlot;
    peano::grid::Grid<exahype::Vertex,exahype::Cell,exahype::State,VertexStack,CellStack,exahype::adapters::PredictionAndPlot2d> _gridWithPredictionAndPlot2d;

     
   exahype::records::RepositoryState               _repositoryState;
   
    tarch::timing::Measurement _measureMeshRefinementCPUTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCPUTime;
    tarch::timing::Measurement _measureInitialConditionAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationCPUTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCPUTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCPUTime;
    tarch::timing::Measurement _measureGridErasingCPUTime;
    tarch::timing::Measurement _measureADERDGTimeStepCPUTime;
    tarch::timing::Measurement _measurePlotAndADERDGTimeStepCPUTime;
    tarch::timing::Measurement _measurePredictionRerunCPUTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCPUTime;
    tarch::timing::Measurement _measureLimiterStatusMergingAndSpreadingMPICPUTime;
    tarch::timing::Measurement _measureLimiterStatusMergingMPICPUTime;
    tarch::timing::Measurement _measureReinitialisationCPUTime;
    tarch::timing::Measurement _measureSolutionRecomputationAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measureNeighbourDataMergingCPUTime;
    tarch::timing::Measurement _measureSolutionUpdateCPUTime;
    tarch::timing::Measurement _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measureTimeStepSizeComputationCPUTime;
    tarch::timing::Measurement _measurePredictionCPUTime;
    tarch::timing::Measurement _measurePredictionAndPlotCPUTime;
    tarch::timing::Measurement _measurePredictionAndPlot2dCPUTime;

    tarch::timing::Measurement _measureMeshRefinementCalendarTime;
    tarch::timing::Measurement _measurePlotAugmentedAMRGridCalendarTime;
    tarch::timing::Measurement _measureInitialConditionAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationCalendarTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationAndPlotCalendarTime;
    tarch::timing::Measurement _measurePredictionAndFusedTimeSteppingInitialisationAndPlot2dCalendarTime;
    tarch::timing::Measurement _measureGridErasingCalendarTime;
    tarch::timing::Measurement _measureADERDGTimeStepCalendarTime;
    tarch::timing::Measurement _measurePlotAndADERDGTimeStepCalendarTime;
    tarch::timing::Measurement _measurePredictionRerunCalendarTime;
    tarch::timing::Measurement _measureLimiterStatusSpreadingCalendarTime;
    tarch::timing::Measurement _measureLimiterStatusMergingAndSpreadingMPICalendarTime;
    tarch::timing::Measurement _measureLimiterStatusMergingMPICalendarTime;
    tarch::timing::Measurement _measureReinitialisationCalendarTime;
    tarch::timing::Measurement _measureSolutionRecomputationAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measureNeighbourDataMergingCalendarTime;
    tarch::timing::Measurement _measureSolutionUpdateCalendarTime;
    tarch::timing::Measurement _measurePostAMRDropMPIMetadataMessagesAndTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measureTimeStepSizeComputationCalendarTime;
    tarch::timing::Measurement _measurePredictionCalendarTime;
    tarch::timing::Measurement _measurePredictionAndPlotCalendarTime;
    tarch::timing::Measurement _measurePredictionAndPlot2dCalendarTime;

   
  public:
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry,
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  computationalDomainOffset
    );
    
    /**
     * Parallel Constructor
     *
     * Used in parallel mode only where the size of the domain is not known 
     * when the type of repository is determined.  
     */
    RepositorySTDStack(
      peano::geometry::Geometry&                   geometry
    );
    
    virtual ~RepositorySTDStack();

    virtual void restart(
      const tarch::la::Vector<DIMENSIONS,double>&  domainSize,
      const tarch::la::Vector<DIMENSIONS,double>&  domainOffset,
      int                                          domainLevel,
      const tarch::la::Vector<DIMENSIONS,int>&     positionOfCentralElementWithRespectToCoarserRemoteLevel
    );
         
    virtual void terminate();
        
    virtual exahype::State& getState();
    virtual const exahype::State& getState() const;
	
    virtual void iterate(int numberOfIterations=1, bool exchangeBoundaryVertices=true);

    virtual void writeCheckpoint(peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> * const checkpoint); 
    virtual void readCheckpoint( peano::grid::Checkpoint<exahype::Vertex, exahype::Cell> const * const checkpoint );
    virtual peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>* createEmptyCheckpoint(); 

    virtual void switchToMeshRefinement();    
    virtual void switchToPlotAugmentedAMRGrid();    
    virtual void switchToInitialConditionAndTimeStepSizeComputation();    
    virtual void switchToPredictionAndFusedTimeSteppingInitialisation();    
    virtual void switchToPredictionAndFusedTimeSteppingInitialisationAndPlot();    
    virtual void switchToPredictionAndFusedTimeSteppingInitialisationAndPlot2d();    
    virtual void switchToGridErasing();    
    virtual void switchToADERDGTimeStep();    
    virtual void switchToPlotAndADERDGTimeStep();    
    virtual void switchToPredictionRerun();    
    virtual void switchToLimiterStatusSpreading();    
    virtual void switchToLimiterStatusMergingAndSpreadingMPI();    
    virtual void switchToLimiterStatusMergingMPI();    
    virtual void switchToReinitialisation();    
    virtual void switchToSolutionRecomputationAndTimeStepSizeComputation();    
    virtual void switchToNeighbourDataMerging();    
    virtual void switchToSolutionUpdate();    
    virtual void switchToPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation();    
    virtual void switchToTimeStepSizeComputation();    
    virtual void switchToPrediction();    
    virtual void switchToPredictionAndPlot();    
    virtual void switchToPredictionAndPlot2d();    

    virtual bool isActiveAdapterMeshRefinement() const;
    virtual bool isActiveAdapterPlotAugmentedAMRGrid() const;
    virtual bool isActiveAdapterInitialConditionAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterPredictionAndFusedTimeSteppingInitialisation() const;
    virtual bool isActiveAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot() const;
    virtual bool isActiveAdapterPredictionAndFusedTimeSteppingInitialisationAndPlot2d() const;
    virtual bool isActiveAdapterGridErasing() const;
    virtual bool isActiveAdapterADERDGTimeStep() const;
    virtual bool isActiveAdapterPlotAndADERDGTimeStep() const;
    virtual bool isActiveAdapterPredictionRerun() const;
    virtual bool isActiveAdapterLimiterStatusSpreading() const;
    virtual bool isActiveAdapterLimiterStatusMergingAndSpreadingMPI() const;
    virtual bool isActiveAdapterLimiterStatusMergingMPI() const;
    virtual bool isActiveAdapterReinitialisation() const;
    virtual bool isActiveAdapterSolutionRecomputationAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterNeighbourDataMerging() const;
    virtual bool isActiveAdapterSolutionUpdate() const;
    virtual bool isActiveAdapterPostAMRDropMPIMetadataMessagesAndTimeStepSizeComputation() const;
    virtual bool isActiveAdapterTimeStepSizeComputation() const;
    virtual bool isActiveAdapterPrediction() const;
    virtual bool isActiveAdapterPredictionAndPlot() const;
    virtual bool isActiveAdapterPredictionAndPlot2d() const;

   
    #ifdef Parallel
    virtual ContinueCommand continueToIterate();
    virtual void runGlobalStep();
    #endif

    virtual void setMaximumMemoryFootprintForTemporaryRegularGrids(double value);
    virtual void logIterationStatistics(bool logAllAdapters) const;
    virtual void clearIterationStatistics();
};


#endif
