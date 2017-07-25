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

#ifndef EXAHYPE_MAPPINGS_Sending_H_
#define EXAHYPE_MAPPINGS_Sending_H_

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
    class Sending;
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
 * <h2>State flags</h2>
 * TODO Currently, this mapping is using the state flag
 * SendMode. We should decompose this mapping into the mappings
 * SendTimeStepData, and SendFaceData
 * in order to minimise the necessity of synchronisation with the global master.
 *
 * @author Dominic Charrier, Tobias Weinzierl
 */
class exahype::mappings::Sending {
public:
  #ifdef Parallel
  /**
   * Allows a rank to skip reduction of time step data if no
   * other data like face unknowns needs to be reduced.
   */
  static bool SkipReductionInBatchedTimeSteps;
  #endif

 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  /**
   * Local copy of the state.
   */
  exahype::State _localState;

  /**
   * A semaphore that is locked if a thread performs a restriction
   * operation.
   */
  static tarch::multicore::BooleanSemaphore _semaphoreForRestriction;

#ifdef Parallel
/**
 * Loop over all the solvers and check
 * if a cell description (ADERDGCellDescription,
 * FiniteVolumesCellDescription,...) is registered
 * for the solver type. If so, send
 * out data or empty messages to the rank \p toRank that
 * owns the neighbouring domain.
 *
 * If not so, send out empty messages for the particular
 * solver.
 *
 * TODO(Dominic): Make messaging solver functionality?
 *
 * \note Not thread-safe.
 */
static void sendSolverDataToNeighbour(
    const int                                    toRank,
    const tarch::la::Vector<DIMENSIONS,int>&     src,
    const tarch::la::Vector<DIMENSIONS,int>&     dest,
    const int                                    srcCellDescriptionIndex,
    const int                                    destCellDescriptionIndex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const int                                    level);

/**
 * Loop over all the solvers and check
 * send out empty messages for the particular
 * solver.
 *
 * \note Not thread-safe.
 */
static void sendEmptySolverDataToNeighbour(
    const int                                     toRank,
    const tarch::la::Vector<DIMENSIONS,int>&      src,
    const tarch::la::Vector<DIMENSIONS,int>&      dest,
    const tarch::la::Vector<DIMENSIONS, double>&  x,
    const int                                     level);
#endif

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
    * Nop.
    */
  Sending();
  #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
  Sending(const Sending& masterThread);
  #endif

  /**
   * If the fine grid cell functions as compute cell for a solver,
   * compute a stable time step size.
   * Further update the time stamp of the compute cell
   * and update the minimum solver time stamp and
   * time step size.
   */
  void enterCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Restricts face data from \p cellDescription to
   * a parent cell description if the fine grid cell associated with
   * cell description is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note We use a semaphore to make this operation thread-safe.
   *
   *
   * <h2> Multicore parallelisation </h2>
   *
   * We face issues here with pfor as the ADER-DG solver can spawn background
   * threads. See the documentation of peano::datatraversal::TasksSet for a
   * remark on this.
   */
  void leaveCell(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
   * Start iteration/grid sweep.
   * Make the state clear its accumulated values.
   *
   * \note Is called once per rank.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Runs over all the registered solvers and sets the
   * reduced minimum time step sizes. Then updates the minimum time stamp
   * of the solvers.
   *
   * \note Is called once per rank.
   */
  void endIteration(exahype::State& solverState);
#ifdef Parallel
  /**
   * Prepare a vertex that is sent to the neighbour
   *
   * We run the $2^d$ adjacent cells and for each cell that is local, we do
   * send its d boundary values that are adjacent to the vertex away.
   *
   * This algorithm translates into a loop very similar to
   * RiemannSolver::touchVertexFirstTime():
   *
   * - Run through all the $2^d$ adjacent cells. Only those that belong to
   *   toRank are of interest. Skip the others. See below for remarks on
   *   interest.
   * - For any cell assigned to toRank, there are d faces that are adjacent to
   *   vertex.
   * - Get the heap indices of all the surrounding cells. Not that some of
   *   them, by definition, are remote.
   *
   * When we run through all the cells adjacent to a vertex, we may communicate
   * only local cells to other ranks. This defines the first two entries in the
   * corresponding if statement. Furthermore, only those cell pairs sharing a
   * face do exchange data. This is done in the third line. The Manhattan
   * distance of the two entries has to be exactly one. Finally (notably on rank
   * 0), we may only send out data if the corresponding cell is inside the
   * domain.
   *
   * <h2>Enumeration</h2>
   *
   * The faces are enumerated: left, right, bottom, top, front, back. Let n be
   * the normal of the face starting with 0. Then, the face index is 2i+f where
   * f means whether it is the face that runs through the left bottom corner (0)
   * or not (1).
   *
   *
   * <h2>MPI administration</h2>
   *
   * Please note that the communication specification deploys the whole heap
   * management to the Peano kernel. There is thus no need to invoke any
   * MPI-specific heap communication operation.
   */
  void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                              const tarch::la::Vector<DIMENSIONS, double>& x,
                              const tarch::la::Vector<DIMENSIONS, double>& h,
                              int level);


  ///////////////////////////////////////
  // WORKER->MASTER Broadcast
  ///////////////////////////////////////
  /**
   * This routine is called on the worker.
   *
   * Send the local array of minimal time step sizes up to the master.
   * Further send up face data if a cell description
   * registered in the fine grid cell of the worker (and master)
   * is of type Ancestor.
   *
   * <h2>Domain Decomposition in Peano</h2>
   * It is important to notice
   * that the master rank's cell
   * and the worker rank's cell
   * overlap.
   *
   * It is further important to notice
   * that both master and worker rank
   * are communicating with their
   * neighbouring cells.
   *
   * The master rank's cell communicates
   * with some neighbour cells
   * via touchVertexFirstTime while
   * the the worker rank's cell communicates
   * with these cells via mergeWithNeighbour,
   * and vice versa.
   *
   * This means that cell descriptions of type
   * Cell registered at a worker cell
   * do not need to communicate their
   * boundary values to their master.
   *
   * Descendants are used to store prolongated
   * face unknowns originating from coarser grid levels.
   * If the overlapping cell holds cell
   * descriptions of type Descendant on the master (and worker) side,
   * we thus need to send data from the master rank to the worker rank.
   * This operation is performed in the mapping Merging since
   * it is a top-down broadcast type operation.
   *
   * Ancestors are used for storing restricted
   * face unknowns originating from finer grid levels.
   * If the overlapping cell holds cell
   * descriptions of type Ancestor on the worker (and master) side,
   * we thus need to send data from the worker rank to the master rank.
   * This operation is performed in the mapping Sending since
   * it is a bottom-up reduction type operation.
   */
  void prepareSendToMaster(
      exahype::Cell& localCell, exahype::Vertex* vertices,
      const peano::grid::VertexEnumerator& verticesEnumerator,
      const exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      const exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

  /**
    * This routine is called on the master.
    *
    * Receive and merge the array of minimal time step sizes send
    * by the worker on the master.
    * Further receive face data if a cell description
    * registered in the fine grid cell is of type Descendant.
    *
    * \see prepareSendToMaster
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
   * Through the result of this routine, we can skip worker-master data
   * transfer as all other mappings return false. Such a skip is advantageous
   * if the runner has decided to trigger multiple grid traversals in one
   * batch. This in turn automatically disables the load balancing.
   *
   * Our strategy thus is as follows: If we may skip the Sending, i.e. the
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
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
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
  void prepareCopyToRemoteNode(exahype::Vertex& localVertex, int toRank,
                               const tarch::la::Vector<DIMENSIONS, double>& x,
                               const tarch::la::Vector<DIMENSIONS, double>& h,
                               int level);

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
   virtual ~Sending();
 #if defined(SharedMemoryParallelisation)
   /**
    * Nop.
    */
   void mergeWithWorkerThread(const Sending& workerThread);
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
