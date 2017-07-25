// This file is part of the Peano project. For conditions of distribution and 
// use, please see the copyright notice at www.peano-framework.org
#ifndef EXAHYPE_ADAPTERS_MeshRefinement2MultiscaleLinkedCell_5_H_
#define EXAHYPE_ADAPTERS_MeshRefinement2MultiscaleLinkedCell_5_H_


#include "tarch/logging/Log.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/MappingSpecification.h"
#include "peano/CommunicationSpecification.h"
#include "peano/grid/VertexEnumerator.h"

#include "exahype/Vertex.h"
#include "exahype/Cell.h"
#include "exahype/State.h"


namespace exahype {
      namespace adapters {
        class MeshRefinement2MultiscaleLinkedCell_5;
      } 
}


/**
 * This is an adapter providing a multiscale linked-cell data structure
 *
 * CellDescriptionsIndex   Name of the index used for the cell indices within the vertex and 
 *          the cell
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.1 $
 */
class exahype::adapters::MeshRefinement2MultiscaleLinkedCell_5 {
  public:
    static peano::MappingSpecification   touchVertexLastTimeSpecification();
    static peano::MappingSpecification   touchVertexFirstTimeSpecification();
    static peano::MappingSpecification   enterCellSpecification();
    static peano::MappingSpecification   leaveCellSpecification();
    static peano::MappingSpecification   ascendSpecification();
    static peano::MappingSpecification   descendSpecification();
    static peano::CommunicationSpecification   communicationSpecification();

    MeshRefinement2MultiscaleLinkedCell_5();

    #if defined(SharedMemoryParallelisation)
    MeshRefinement2MultiscaleLinkedCell_5(const MeshRefinement2MultiscaleLinkedCell_5& masterThread);
    #endif

    virtual ~MeshRefinement2MultiscaleLinkedCell_5();
  
    #if defined(SharedMemoryParallelisation)
    void mergeWithWorkerThread(const MeshRefinement2MultiscaleLinkedCell_5& workerThread);
    #endif

    void createInnerVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void createBoundaryVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    /**
     * @todo Only half of the documentation
     * @todo Enumeration has changed
     *
     * In an adaptive grid, not all of the $2^d$ adjacent cells exist for hanging
     * vertices. Since each vertex is supposed to hold the adjacent vertices in
     * order to fill the ghostlayers of the cellDescriptiones appropriately, the adjacent
     * indices of hanging vertices need to be filled by the data of the vertices
     * on the next coarser grid. This filling is implemented in this method.
     *
     * !!! The idea
     * Each vertex holds $2^d$ indices. In the vertices they are numbered from 0
     * to $2^d-1$. However, in this method they are considered to exist in a
     * n-dimensional array. In 2d this would look like
     *
     * (0,1)|(1,1)
     * -----v-----
     * (0,0)|(1,0)
     *
     * The linearization looks as follow:
     *
     *   1  |  0
     * -----v-----
     *   3  |  2
     *
     * In the following the term "fine grid" refers to the $4^d$ vertices
     * belonging to the $3^d$ fine grid cells which overlap with the coars grid
     * cell.
     *
     * On the coarse grid cell we again consider the vertices being arranged in a
     * n-dimensional array:
     *
     * (0,1)-----(1,1)
     *   |          |
     *   |          |
     *   |          |
     * (0,0)-----(1,0)
     *
     * Each of them hold again the $2^d$ adjacent indices, while those which refer
     * to a refined cell are set to -1. A hanging vertex therefore gets the
     * adjacent indices from the nearest coarse grid vertex. If they coincide the
     * data can just be used directly. If not, it depends on which boundary of the
     * coarse grid cell the hanging vertex resides. Here the (single) index
     * outside of the coarse grid cell is assigned for all indices of the hanging
     * vertex pointing in the direction of this neighboring coarse grid cell.
     *
     * !!! The algorithm
     * It gets a hanging vertex and performs a loop over the $2^d$ adjacent-cellDescription-
     * indices.
     * In each loop iteration it computes the n-dimensional index of the coarse
     * grid vertex (fromCoarseGridVertex) from which the data has to be copied.
     * For each dimension d with $0\le d <n$:
     *  - If the fine grid position of the hanging vertex in dimension $d$ is 0 set
     *    $fromCoarseGridVertex(d)$ to 0. If it is equals 3 then set
     *    $fromCoarseGridVertex(d)$ to 1. By this we ensure that we always choose
     *    the nearest coarse grid vertex in dimension $d$. If the hanging vertex
     *    resides in a corner of the fine grid this approach always chooses the
     *    coarse grid vertex that is located on the same position.
     *  - If the fine grid position of the hanging vertex in dimension $d$ is
     *    neither 0 nor 3 then the value of $fromCoarseGridVertex(d)$ depends on
     *    the adjacent-cellDescription-index $k$ that has to be set currently. $k(d)$ can
     *    either be 0 or 1. If $k(d)$ is 0 than we want to get data from the
     *    in this dimension "lower" coarse grid vertex, so we set
     *    $fromCoarseGridVertex(d)$ to 0 as well. In the case of $k(d)=1$ we set
     *    $fromCoarseGridVertex(d)$ to 1, accordingly. This actually doesn't
     *    matter since the appropriate adjacent-cellDescription-indices of the to coarse
     *    grid vertices have to be the same, since they are pointing to the same
     *    adjacent cell.
     * The determination of the correct adjacent-cellDescription-index of the coarse grid
     * vertex (coarseGridVertexAdjacentCellDescriptionIndex) is done in a similar way. So,
     * for the adjacent-cellDescription-index $k$ on the hanging vertex:
     *  - As stated before, if the fine and coarse grid vertices coincide we can
     *    just copy the adjacent-cellDescription-index. Therefore, if the fine grid position
     *    of the hanging vertex in dimension $d$ is equal to 0 or to 3, we set
     *    $coarseGridVertexAdjacentCellDescriptionIndex(d)$ to $k(d)$.
     *  - Otherwise, we just set $coarseGridVertexAdjacentCellDescriptionIndex(d)$ to the
     *    inverted $k(d)$. I.e. if $k(d) = 0$ we set
     *    $coarseGridVertexAdjacentCellDescriptionIndex(d)$ to 1 and the other way around.
     *
     */
    void createHangingVertex(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void destroyHangingVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void destroyVertex(
      const exahype::Vertex&   fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );


    void createCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const         fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void destroyCell(
      const exahype::Cell&           fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );
        
    #ifdef Parallel
    void mergeWithNeighbour(
      exahype::Vertex&  vertex,
      const exahype::Vertex&  neighbour,
      int                                           fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareSendToNeighbour(
      exahype::Vertex&  vertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      exahype::Vertex&  localVertex,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );

    void prepareCopyToRemoteNode(
      exahype::Cell&  localCell,
      int                                           toRank,
      const tarch::la::Vector<DIMENSIONS,double>&   cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&   cellSize,
      int                                           level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Vertex&  localVertex,
      const exahype::Vertex&  masterOrWorkerVertex,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  x,
      const tarch::la::Vector<DIMENSIONS,double>&  h,
      int                                       level
    );

    void mergeWithRemoteDataDueToForkOrJoin(
      exahype::Cell&  localCell,
      const exahype::Cell&  masterOrWorkerCell,
      int                                       fromRank,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                       level
    );

    bool prepareSendToWorker(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell,
      int                                                                  worker
    );

    void prepareSendToMaster(
      exahype::Cell&                       localCell,
      exahype::Vertex *                    vertices,
      const peano::grid::VertexEnumerator&       verticesEnumerator, 
      const exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&       coarseGridVerticesEnumerator,
      const exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&   fineGridPositionOfCell
    );

    void mergeWithMaster(
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
    );


    void receiveDataFromMaster(
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
    );


    void mergeWithWorker(
      exahype::Cell&           localCell, 
      const exahype::Cell&     receivedMasterCell,
      const tarch::la::Vector<DIMENSIONS,double>&  cellCentre,
      const tarch::la::Vector<DIMENSIONS,double>&  cellSize,
      int                                          level
    );


    void mergeWithWorker(
      exahype::Vertex&        localVertex,
      const exahype::Vertex&  receivedMasterVertex,
      const tarch::la::Vector<DIMENSIONS,double>&   x,
      const tarch::la::Vector<DIMENSIONS,double>&   h,
      int                                           level
    );
    #endif


    void touchVertexFirstTime(
      exahype::Vertex&               fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                          fineGridH,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfVertex
    );


    void touchVertexLastTime(
      exahype::Vertex&         fineGridVertex,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridX,
      const tarch::la::Vector<DIMENSIONS,double>&                    fineGridH,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                       fineGridPositionOfVertex
    );
    

    void enterCell(
      exahype::Cell&                 fineGridCell,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&                             fineGridPositionOfCell
    );


    void leaveCell(
      exahype::Cell&                          fineGridCell,
      exahype::Vertex * const                 fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const                 coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&                          coarseGridCell,
      const tarch::la::Vector<DIMENSIONS,int>&      fineGridPositionOfCell
    );


    void beginIteration(
      exahype::State&  solverState
    );


    void endIteration(
      exahype::State&  solverState
    );

    void descend(
      exahype::Cell * const          fineGridCells,
      exahype::Vertex * const        fineGridVertices,
      const peano::grid::VertexEnumerator&                fineGridVerticesEnumerator,
      exahype::Vertex * const        coarseGridVertices,
      const peano::grid::VertexEnumerator&                coarseGridVerticesEnumerator,
      exahype::Cell&                 coarseGridCell
    );


    void ascend(
      exahype::Cell * const    fineGridCells,
      exahype::Vertex * const  fineGridVertices,
      const peano::grid::VertexEnumerator&          fineGridVerticesEnumerator,
      exahype::Vertex * const  coarseGridVertices,
      const peano::grid::VertexEnumerator&          coarseGridVerticesEnumerator,
      exahype::Cell&           coarseGridCell
    );    
};


#endif
