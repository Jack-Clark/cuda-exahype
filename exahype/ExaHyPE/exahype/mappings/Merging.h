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

#ifndef EXAHYPE_MAPPINGS_Merging_H_
#define EXAHYPE_MAPPINGS_Merging_H_

#include "tarch/la/Vector.h"
#include "tarch/logging/Log.h"

#include "peano/CommunicationSpecification.h"
#include "peano/MappingSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "tarch/multicore/MulticoreDefinitions.h"

#include "exahype/Cell.h"
#include "exahype/State.h"
#include "exahype/Vertex.h"

namespace exahype {
  namespace mappings {
    class Merging;
  }
}

/**
 * Run the Riemann solves on the faces
 *
 * @todo Dominic, bitte ordentlich dokumentieren, was hier wann, wo und warum passiert.
 * Bitte auch dokumentieren, falls Du was mal probiert hast und es nicht funktioniert hat.
 *
 * <h2>State flags</h2>
 * TODO Currently, this mapping is using the state flag MergeMode
 * We should decompose this mapping into the mappings MergeTimeStepData,
 * and MergeTimeStepAndFaceData in order to minimise the necessity of
 * synchronisation with the global master.
 *
 * <h2>Synchronisation in multi-rank environment</h2>
 *
 * The individual Riemann solves have to have knowledge about
 * admissible time step sizes. This information is read from the solver
 * instances. In return, they update the solver instances with the maximal
 * permitted time step size for the subsequent iteration. This synchronisation
 * is realised in the mapping NewTimeStep, i.e. no solver synchronisation is to
 * be done in this mapping. It is also realised in SpaceTimePredictor. So
 * please consult the documentation there.
 *
 * This mapping is used to synchronise the compute cell time step sizes
 * with the solver ones. It further is used to reset the Riemann
 * solve flags on compute and helper cells.
 *
 * @author Dominic E. Charrier and Tobias Weinzierl
 */
class exahype::mappings::Merging {
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
   * Temporary variable per solver for storing
   * space-time face unknowns.
   */
//  double**  _tempSpaceTimeFaceUnknownsArray  = nullptr; todo

  /**
   * Temporary variable per solver for storing
   * face unknowns.
   */
  double***  _tempFaceUnknowns = nullptr;

  /**
   * Temporary variables per solver for storing state sized (=number of variables)
   * quantities like eigenvalues or averaged states.
   */
  double*** _tempStateSizedVectors = nullptr;

  /**
   * Temporary variable per solver for storing square matrices
   * of the size number of variables times number of variables.
   */
  double*** _tempStateSizedSquareMatrices = nullptr;

  /**
   * Initialises temporary variables.
   */
  void prepareTemporaryVariables();

  /**
   * Free memory allocated for temporary variables.
   */
  void deleteTemporaryVariables();

  #ifdef Debug // TODO(Dominic): Exclude shared memory etc.
  /*
   *  Counter for the interior face solves for debugging purposes.
   */
  int _interiorFaceMerges;
  /*
   *  Counter for the boundary face solves for debugging purposes.
   */
  int _boundaryFaceMerges;
  #endif

  #ifdef Parallel
  /**
   * Iterates over the received metadata and every time
   * we find a valid entry, we call mergeWithNeighbourData
   * on the solver corresponding to the metadata.
   * if we want to receive the neighbour data
   * or if we just want to drop it.
   * \note Not thread-safe.
   */
  void mergeWithNeighbourData(
      const int fromRank,
      const int srcCellDescriptionIndex,
      const int destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);

  /**
   * Iterates over the received metadata and
   * drops the received neighbour data.
   *
   * \note Not thread-safe.
   */
  static void dropNeighbourData(
      const int fromRank,
      const int srcCellDescriptionIndex,
      const int destCellDescriptionIndex,
      const tarch::la::Vector<DIMENSIONS,int>& src,
      const tarch::la::Vector<DIMENSIONS,int>& dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int level,
      const exahype::MetadataHeap::HeapEntries& receivedMetadata);
  #endif

public:

  /**
   * Call the touch vertex first time event on every vertex of
   * the grid. Run in parallel but avoid fine grid races.
   *
   * <h2>Shared Memory</h2>
   * The AvoidFineGridRaces multithreading specification prevents that
   * more than one threads write data for the same face of the grid
   * at the same time.
   *
   * The specification is realised by touching the vertices in
   * a red-black (X and O in the figure below) manner:
   * We might first process the X-vertices in parallel
   * and then the O-vertices.
   * In each step, faces adjacent to the vertices,
   * do not overlap and race conditions can thus
   * not occur.
   *
   *     |    |
   *     |    |
   * ----O----X-----
   *     |    |
   *     |    |
   * ----X----O-----
   *     |    |
   *     |    |
   *
   * TODO(Dominic): It might be useful to introduce a multithreading specification
   * "AvoidFineGridRacesOnlyRed" that processes only the red
   * vertices and not the black ones. In fact, the second sweep over the black vertices only
   * finds the riemannSolvePerfomed flags set and does nothing in
   * our current implementation.
   */
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  /**
   * Nop.
   */
  static peano::MappingSpecification enterCellSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  /**
   * Nop
   */
  static peano::MappingSpecification leaveCellSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification ascendSpecification();

  /**
   * Nop
   */
  static peano::MappingSpecification descendSpecification();

  /**.
   * The mapping does synchronise through synchroniseTimeStepping() invoked
   * on the solvers. Yet, no data is transported through the vertices, the
   * cell or the state object.
   *
   * Though we need valid solvers for the actual Riemann solves, we do not
   * do any solver exchange in this mapping. The appropriate data exchange
   * is done in NewTimeStep() and SolutionUpdate().
   *
   * Further let Peano handle heap data exchange internally.
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
   * Solve Riemann problems on all interior faces that are adjacent
   * to this vertex and impose boundary conditions on faces that
   * belong to the boundary. This is done for all cell descriptions
   * belonging to the cells that are an interior face.
   *
   * The routine itself runs the loop over the faces. The actual
   * functionality is outsourced to solveRiemannProblemAtInterface().
   *
   * The function ensures implicitly that interior faces
   * do not align with MPI boundaries. In this case, no operation
   * is performed.
   *
   * This method sets the riemannSolvePerformed flag on a cell description
   * if boundary conditions have been imposed for this cell description.
   * This method sets the riemannSolvePerformed flag on both cell descriptions
   * (per solver) for interior faces if a Riemann solve has been performed for
   * both cell descriptions.
   *
   * \note The function itself is not thread-safe.
   * Thread-safety of this function is ensured by setting
   * RiemannSolver::touchVertexFirstTimeSpecification()
   * to peano::MappingSpecification::AvoidFineGridRaces.
   *
   * <h2>Limiter identification</h2>
   * Each ADER-DG solver analyses the local min and max values within a cell.
   * This information however is not stored in the cell but on the 2d faces
   * of a cell. See Cell::setSolutionMinMaxAndAnalyseValidity() for details.
   * The face information then in the subsequent step has to be merged which
   * is done when we trigger the Riemann solve.
   *
   * @see Cell::mergeSolutionMinMaxOnFace()
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
   * TODO(Dominic): Add docu.
   */
  Merging();

#if defined(SharedMemoryParallelisation)
  /**
   * TODO(Dominic): Add docu.
   */
  Merging(const Merging& masterThread);
#endif

  /**
   * Free previously allocated temporary variables.
   */
  virtual ~Merging();

  /**
   * Initialise temporary variables
   * if they are not initialised yet (or
   * if a new solver was introuced to the grid.
   * This is why we put the initialisation
   * in beginIteration().
   *
   * In debug mode, further resets counters for Riemann solves at
   * interior and boundary faces.
   */
  void beginIteration(exahype::State& solverState);

  /**
   * Frees previously allocated temporary variables.
   *
   * In debug mode, prints the output of counters.
   */
  void endIteration(exahype::State& solverState);

#if defined(SharedMemoryParallelisation)
  /**
   * Nop.
   */
  void mergeWithWorkerThread(const Merging& workerThread);
#endif


#ifdef Parallel
  /**
   * Merge vertex with the incoming vertex from a neighbouring computation node.
   *
   * When Peano is running in parallel the data exchange is done vertex-wise
   * between two grid iterations, i.e. the predictor sends out data in one step
   * and in the following step we receive this data and merge it into the local
   * Riemann integrals. The routine is thus very similar to
   * touchVertexFirstTime():
   *
   * - We identify incoming faces.
   * - We create (temporary) indices on the heap.
   * - We receive the data into these indices.
   *
   * As we always receive data in the iteration following the sends and as
   * Peano inverts the traversal direction after each grid sweep, we have to
   * invert the order in which data is received, too.
   *
   * \note It happens or is possible that this operation is performed
   * after touchVertexLastTime(...) was invoked.
   *
   * <h2>Send</h2>
   * Sending out data corresponds logically to a projection of cell data onto
   * the faces. Therefore, I realise it within
   * SpaceTimePredictor::prepareSendToNeighbour().
   *
   *
   * <h2>Min max analysis</h2>
   * The min/max analysis runs analogously. We do send out min and max from
   * either side to the other rank and then merge min and max on both sides
   * into the local data.
   *
   * <h2>Face data exchange counters</h2>
   * Cell::isInside() does not imply that all adjacent vertices are
   * inside. If we count down the counter only on
   * vertices that are inside we might not send out all faces
   * of a cell that is close to the boundary.
   * @see Prediction::countListingsOfRemoteRankByInsideVerticesAtFace.
   */
  void mergeWithNeighbour(exahype::Vertex& vertex,
                          const exahype::Vertex& neighbour, int fromRank,
                          const tarch::la::Vector<DIMENSIONS, double>& x,
                          const tarch::la::Vector<DIMENSIONS, double>& h,
                          int level);

  /**
   * This routine is called on the master.
   *
   * Send the local array of reduced minimal time step sizes down to
   * the worker. Further send down face data if a cell description
   * registered in the fine grid cell is of type Descendant.
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
  bool prepareSendToWorker(
      exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      int worker);

  /**
   * Receive kick-off message from master
   *
   * Counterpart of prepareSendToWorker(). This operation is called once when
   * we receive data from the master node.
   *
   * To do a proper Riemann solve, it is important that all the solvers have
   * the right state.We therefore receive for each individual solver a couple
   * of messages and merge them into the global instance before we continue
   * with the actual iteration. However, this data exchange in the mapping
   * Sending. See the documentation there for rationale.
   *
   * We further receive face data send from the master in
   * this hook if a cell description of the worker (and master)
   * is of type Descendant. See prepareSendToWorker(...) for a rationale
   *
   * \see prepareSendToWorker(...)
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

  //
  // Below all methods are nop.
  //
  //===================================


  /**
   * Nop.
   */
  void prepareSendToNeighbour(
      exahype::Vertex& vertex, int toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const tarch::la::Vector<DIMENSIONS, double>& h, int level);

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
