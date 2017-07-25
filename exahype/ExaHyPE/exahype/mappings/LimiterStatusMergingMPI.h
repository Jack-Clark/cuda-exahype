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

#ifndef EXAHYPE_MAPPINGS_LimiterStatusMergingMPI_H_
#define EXAHYPE_MAPPINGS_LimiterStatusMergingMPI_H_

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
class LimiterStatusMergingMPI;
}
}

/**
 * This mapping is part of three mappings (+adapters) which together perform the limiter
 * status spreading the recomputation of troubled cells (and their direct) neighbours.
 *
 * \see exahype::mappings::SolutionRecomputation for more details.
 *
 * @author Dominic Charrier
 */
class exahype::mappings::LimiterStatusMergingMPI {
 private:
  /**
   * Logging device for the trace macros.
   */
  static tarch::logging::Log _log;

  #ifdef Debug
  int _interiorFaceMerges;
  int _boundaryFaceMerges;
  #endif

 public:
  /**
   * Run through the whole grid. Run concurrently on the fine grid.
   */
  static peano::MappingSpecification enterCellSpecification();

  /**
   * Nop.
   */
  static peano::MappingSpecification touchVertexLastTimeSpecification();
  /**
   * Nop.
   */
  static peano::MappingSpecification touchVertexFirstTimeSpecification();
  /**
   * Nop.
   */
  static peano::MappingSpecification leaveCellSpecification();
  /**
   * Nop.
   */
  static peano::MappingSpecification ascendSpecification();
  /**
   * Nop.
   */
  static peano::MappingSpecification descendSpecification();

  /**
   * Mask out data exchange between master and worker.
   * Further let Peano handle heap data exchange internally.
   */
  static peano::CommunicationSpecification communicationSpecification();

  /**
     * Initialise debug counters.
     */
    void beginIteration(exahype::State& solverState);

    /**
     * Print debug counters.
     */
    void endIteration(exahype::State& solverState);

    /**
     * Overwrite the merged limiter status values on the faces
     * by a unified value determined on hand of those values.
     * Do not touch the cell-based limiter status.
     */
    void enterCell(
        exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell);

    /**
     * TODO(Dominic): Add docu.
     */
    void touchVertexFirstTime(
        exahype::Vertex& fineGridVertex,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
        const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
        exahype::Vertex* const coarseGridVertices,
        const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
        exahype::Cell& coarseGridCell,
        const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex);

  #ifdef Parallel
    /**
     * TODO(Dominic): Add docu.
     */
    void mergeWithNeighbour(exahype::Vertex& vertex,
                            const exahype::Vertex& neighbour, int fromRank,
                            const tarch::la::Vector<DIMENSIONS, double>& x,
                            const tarch::la::Vector<DIMENSIONS, double>& h,
                            int level);

    /**
     * We only drop the received limiter status for LimitingADERDGSolvers
     * where we have detected a change of the limiter domain.
     * This information should be available on all ranks.
     * We ignore other solver types.
     *
     * TODO(Dominic): Add more docu.
     */
    static void dropNeighbourMergedLimiterStatus(
        const int                                    fromRank,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const int                                    srcCellDescriptionIndex,
        const int                                    destCellDescriptionIndex,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level,
        const exahype::MetadataHeap::HeapEntries&    receivedMetadata);

    /**
     * We only merge the face-wise limiter status for LimitingADERDGSolvers
     * where we have detected a change of the limiter domain.
     * This information should be available on all ranks.
     * We ignore other solver types.
     *
     * TODO(Dominic): Add more docu.
     */
    static void mergeNeighourMergedLimiterStatus(
        const int                                    fromRank,
        const tarch::la::Vector<DIMENSIONS,int>&     src,
        const tarch::la::Vector<DIMENSIONS,int>&     dest,
        const int                                    srcCellDescriptionIndex,
        const int                                    destCellDescriptionIndex,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level,
        const exahype::MetadataHeap::HeapEntries&    receivedMetadata);

    /**
     * TODO(Dominic): Add docu.
     */
    void prepareSendToNeighbour(exahype::Vertex& vertex, int toRank,
                                const tarch::la::Vector<DIMENSIONS, double>& x,
                                const tarch::la::Vector<DIMENSIONS, double>& h,
                                int level);



    //
    // Below all methods are nop.
    //
    //===================================

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
     * Nop
     */
    LimiterStatusMergingMPI();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    LimiterStatusMergingMPI(const LimiterStatusMergingMPI& masterThread);
  #endif

    /**
     * Nop.
     */
    virtual ~LimiterStatusMergingMPI();

  #if defined(SharedMemoryParallelisation)
    /**
     * Nop.
     */
    void mergeWithWorkerThread(const LimiterStatusMergingMPI& workerThread);
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
