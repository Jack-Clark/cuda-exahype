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
 
#ifndef EXAHYPE_MAPPINGS_FusedTimeSteppingInitialisation_H_
#define EXAHYPE_MAPPINGS_FusedTimeSteppingInitialisation_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/BooleanSemaphore.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

namespace exahype {
  namespace mappings {
    class FusedTimeSteppingInitialisation;
  }
}

/**
 * This mapping is used to set the corrector time stamp and
 * step size to the predictor values for every ADER-DG solver
 * itself and the cell descriptions belonging to this solver.
 * Further adds the predictor time step
 * size on the predictor time stamp for
 * every ADER-DG solver and its
 * cell descriptions.
 *
 * !!! Rationale:
 * In the fused time stepping scheme, the predictor
 * time stamp must be one time step size ahead of the
 * corrector time stamp.
 *
 * @author Dominic Charrier
 */
class exahype::mappings::FusedTimeSteppingInitialisation {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * Sets the corrector time stamp and step size to the
   * predictor values. Adds the predictor time step
   * size on the predictor time stamp.
   *
   * !!! Rationale:
   * In the fused time stepping scheme, the predictor
   * time stamp must be one time step size ahead of the
   * corrector time stamp.
   */
  void initialiseFusedTimestepping(exahype::solvers::Solver* solver) const;

  /**
   * Similar to ::initialiseFusedTimestepping(exahype::solvers::Solver*) but per cell description.
   */
  void initialiseFusedTimestepping(exahype::solvers::Solver* solver,const int cellDescriptionsIndex, const int element) const;
 public:
  /**
   * Run through whole tree. Run concurrently on fine grid.
   */
  static peano::MappingSpecification enterCellSpecification();

  /**
   * Nop.
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  static peano::MappingSpecification leaveCellSpecification();
  static peano::MappingSpecification ascendSpecification();
  static peano::MappingSpecification descendSpecification();

  /**
   * TODO(Dominic): Currently, we need to broadcast the state in order
   * to get the information if the fused time stepping algorithm
   * was chosen by the user. However since we agreed that
   * status does not change during the simulation time,
   * we can set a on each rank's State a static flag
   * in Runner::run();
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
   * \see ::initialiseFusedTimestepping(exahype::solvers::Solver*,const int,const int).
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  virtual ~FusedTimeSteppingInitialisation();
  #if defined(SharedMemoryParallelisation)
  /**
   * Copy the state from the master.
   */
  FusedTimeSteppingInitialisation(const FusedTimeSteppingInitialisation& masterThread);
  #endif

  /**
   * \see ::initialiseFusedTimestepping(exahype::solvers::Solver*).
   */
  void beginIteration(exahype::State& solverState);


  //
  // Below every routine is nop.
  //
  // ==================================


  /**
   * Nop.
   */
  void endIteration(exahype::State& solverState);
#ifdef Parallel
  /**
   * Nop.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Nop.
   */
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);

  /**
   * Nop.
   */
  void mergeWithMaster(
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
      exahype::State& masterState);

  /**
   * Nop.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);
  /**
   * Nop.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);
  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);
  /**
   * Nop.
   */
  void prepareCopyToRemoteNode(
      exahype::Cell& localCell, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);
  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);
  /**
   * Nop.
   */
  void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
      int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level);
  /**
   * Nop.
   */
  void receiveDataFromMaster(
      exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
      const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
      exahype::Vertex* const receivedCoarseGridVertices,
      const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
      exahype::Cell& receivedCoarseGridCell,
      exahype::Vertex* const workersCoarseGridVertices,
      const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
      exahype::Cell& workersCoarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Cell& localCell,
                       const exahype::Cell& receivedMasterCell,
                       const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                       const tarch::la::Vector<DIMENSIONS, double>& cellSize,
                       int level);
  /**
   * Nop.
   */
  void mergeWithWorker(exahype::Vertex& localVertex,
                       const exahype::Vertex& receivedMasterVertex,
                       const tarch::la::Vector<DIMENSIONS, double>& x,
                       const tarch::la::Vector<DIMENSIONS, double>& h,
                       int level);
#endif
  /**
   * Nop.
   */
  FusedTimeSteppingInitialisation();

 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   void mergeWithWorkerThread(const FusedTimeSteppingInitialisation& workerThread);
 #endif
   /**
    * Nop.
    */
   void createInnerVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createBoundaryVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createHangingVertex(
       exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void destroyHangingVertex(
       const exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void destroyVertex(
       const exahype::Vertex& fineGridVertex,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
       const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
   /**
    * Nop.
    */
   void createCell(
       exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
       const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
   /**
    * Nop.
    */
   void destroyCell(
       const exahype::Cell& fineGridCell,
       exahype::Vertex* const fineGridVertices,
       const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
       exahype::Vertex* const coarseGridVertices,
       const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
       exahype::Cell& coarseGridCell,
       const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void touchVertexFirstTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void touchVertexLastTime(
      exahype::Vertex& fineGridVertex,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
      const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);
  /**
   * Nop.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);
  /**
   * Nop.
   */
  void descend(
      exahype::Cell* const fineGridCells,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell);
  /**
   * Nop.
   */
  void ascend(exahype::Cell* const fineGridCells,
              exahype::Vertex* const fineGridVertices,
              const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
              exahype::Vertex* const coarseGridVertices,
              const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
              exahype::Cell& coarseGridCell);
};
#endif
