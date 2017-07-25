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
 
#ifndef EXAHYPE_MAPPINGS_TimeStepSizeComputation_H_
#define EXAHYPE_MAPPINGS_TimeStepSizeComputation_H_

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
    class TimeStepSizeComputation;
  }
}

/**
 * Determine a global time step size
 *
 * The global time step computation runs through all the cells. Per cell, it
 * runs through all involved solvers and determines the corresponding minimal
 * time step sizes. Once the traversal terminates, all solvers thus know what
 * the minimal permitted time step size is. We now take the minimal time step
 * sizes and inform the solvers about them through
 * updateMinNextPredictorTimeStepSize().
 *
 * In the subsequent time step, these minimal time step sizes then are used by
 * synchroniseTimeStepping() (see notably the mapping NewTimeStep) to move the
 * patch forward in time. There is no need to take extra care for the
 * optimistic time step choice - we can determine from outside whether we tend
 * to overshoot and thus have to rerun the predictor. This is done in the
 * runner.
 *
 * <h2>Multicore parallelisation</h2>
 * See documentation of _minTimeStepSizes/
 *
 * <h2>MPI parallelisation</h2>
 *
 * @author Dominic Charrier, Tobias Weinzierl
 */
class exahype::mappings::TimeStepSizeComputation {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * A minimum time step size for each solver.
   */
  std::vector<double> _minTimeStepSizes;

  /**
   * A minimum cell size for each solver.
   */
  std::vector<double> _minCellSizes;

  /**
   * A maximum cell size for each solver.
   */
  std::vector<double> _maxCellSizes;

  double** _tempEigenValues = nullptr;

  /**
   * Prepare a appropriately sized vector _minTimeStepSizes
   * with elements initiliased to MAX_DOUBLE.
   */
  void prepareLocalTimeStepVariables();

  /**
   * Per solver, allocate a temporary eigenvalues
   * array.
   */
  void prepareTemporaryVariables();

  /**
   * Free memory reserved for the eigenvalue vectors we have
   * allocated per solver.
   */
  void deleteTemporaryVariables();

  /**
   * Reinitialises the corrector and predictor time step sizes of an ADER-DG solver
   * with stable time step sizes if we detect a-posteriori that the CFL condition was
   * harmed by the estimated predictor time step size used in the last iteration.
   */
  void reinitialiseTimeStepDataIfLastPredictorTimeStepSizeWasInstable(exahype::State& state,exahype::solvers::Solver* solver) const;

  /**
   * If the original time stepping algorithm is used for the ADER-DG scheme,
   * we need to enforce that the corrector time step size is identical to the
   * predictor time step size.
   * We further need to
   */
  void reconstructStandardTimeSteppingData(exahype::solvers::Solver* solver) const;


  /**
   * Similar to ::overwriteCorrectorTimeStepDataWithPredictorValues(exahype::solvers::Solver*) but for
   * a cell description.
   */
  void reconstructStandardTimeSteppingData(exahype::solvers::Solver* solver,const int cellDescriptionsIndex,const int element) const;

 public:
  /**
   * This flag is used by mapping InitialCondition to veto the
   * fused time step size reinitialisation after the computation
   * of the first time step size.
   *
   * TODO(Dominic): Is that okay? Or is it a hack?
   */
  static bool VetoFusedTimeSteppingTimeStepSizeReinitialisation;

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
   * The global time step computation does synchronise the individual cells
   * with the solver instances from the master. The actual
   * synchronisation/consistency routines
   * for the solvers are done in NewTimeStep and
   * the SpaceTimePredictor.
   *
   * The fundamental job of the global time step mapping is to report back to
   * the master what time step is permitted. As all those operations are
   * actually done in enterCell---also the veto of a global time step is done
   * in the Riemann solver, i.e. before enterCell---we can send back data as
   * soon as the traversal operation leaves the local subtree.
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
   * If the fine grid cell functions as compute cell for a solver,
   * compute a stable time step size.
   *
   * Then update the time stamp of the compute cell
   * and update the minimum solver time stamp and
   * time step size.
   *
   * Finally, update the minimum and maximum mesh cell size
   * fields per solver with the size of the fine grid cell.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Delete temporary variables for all threads.
   */
  virtual ~TimeStepSizeComputation();
  #if defined(SharedMemoryParallelisation)
  /**
   * Initialise temporary variables for worker threads.
   */
  TimeStepSizeComputation(const TimeStepSizeComputation& masterThread);
  #endif

  /**
   * Start iteration/grid sweep.
   * Make the state clear its accumulated values.
   *
   * Further initialise temporary variables
   * if they are not initialised yet (or
   * if a new solver was introuced to the grid.
   * This is why we put the initialisation
   * in beginIteration().
   *
   * \note Is called once per rank.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Runs over all the registered solvers and sets the
   * reduced minimum time step sizes. Then updates the minimum time stamp
   * of the solvers.
   *
   * Iterate over the solvers and start a new time step
   * on every solver.
   *
   * <h2>Fused ADER-DG time stepping</h2>
   * If we use the fused ADER-DG time stepping algorithm,
   * The solver or (the solver belonging to the global master in the MPI context)
   * is not allowed to perform the time step update
   * directly. It first has to check if the previously used
   * min predictor time step size was stable one.
   * Otherwise, we would corrupt the corrector time stamp
   * with an invalid value.
   *
   * <h2>MPI</h2>
   * Here we start again a new time step "in the small" on the
   * worker rank and overwrite it later on again if a synchronisation is applied
   * by the master rank.
   *
   * It is important to keep in mind that endIteration() on a worker
   * is called before the prepareSendToMaster routine.
   * We thus send out the current time step size from
   * the worker to the master.
   *
   * On the master, the mergeWithMaster routine is however called
   * before endIteration.
   * We thus merge the received time step size with the next
   * time step size on the master.
   *
   * \see exahype::mappings::Sending,exahype::mappings::Merging
   */
  void endIteration(exahype::State& solverState);
#ifdef Parallel
  /**
   * This routine is called on the worker.
   *
   * Send the local array of minimal time step sizes up to the master. This is
   * one MPI_Send on the whole array.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Through the result of this routine, we can skip worker-master data
   * transfer as all other mappings return false. Such a skip is advantageous
   * if the runner has decided to trigger multiple grid traversals in one
   * batch. This in turn automatically disables the load balancing.
   *
   * Our strategy thus is as follows: If we may skip the reduction, i.e. the
   * user has enabled this optimisation in the ExaHyPE spec file, then we
   * return false if load balancing is disabled.
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
   * This routine is called on the master.
   * TODO(Dominic): Docu.
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
  TimeStepSizeComputation();

 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   void mergeWithWorkerThread(const TimeStepSizeComputation& workerThread);
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
