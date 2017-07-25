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
 
#ifndef _EXAHYPE_CELL_H_
#define _EXAHYPE_CELL_H_

#include "exahype/State.h"

#include "exahype/records/Cell.h"
#include "peano/grid/Cell.h"
#include "peano/grid/VertexEnumerator.h"

#include "exahype/records/ADERDGCellDescription.h"
#include "exahype/records/FiniteVolumesCellDescription.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

namespace exahype {
  class Cell;

  namespace solvers {
    class ADERDGSolver;
    class FiniteVolumesSolver;
  }
}

/**
 * @todo Dominic We should add some proper descriptions here one day.
 */
class exahype::Cell : public peano::grid::Cell<exahype::records::Cell> {
 private:
  typedef class peano::grid::Cell<exahype::records::Cell> Base;

  static tarch::logging::Log _log;

 public:
  /**
   * Default Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it.
   */
  Cell();

  /**
   * This constructor should not set any attributes. It is used by the
   * traversal algorithm whenever it allocates an array whose elements
   * will be overwritten later anyway.
   */
  Cell(const Base::DoNotCallStandardConstructor&);

  /**
   * Constructor
   *
   * This constructor is required by the framework's data container. Do not
   * remove it. It is kind of a copy constructor that converts an object which
   * comprises solely persistent attributes into a full attribute. This very
   * functionality is implemented within the super type, i.e. this constructor
   * has to invoke the correponsing super type's constructor and not the super
   * type standard constructor.
   */
  Cell(const Base::PersistentCell& argument);

  /**
   * Returns if the face is inside. Inside for
   * a face means that at least one vertex
   * of the 2^{d-1} vertices building up
   * the face must be inside of the domain.
   *
   * Otherwise the face might be on the
   * boundary of the domain or outside of
   * the domain.
   */
  static bool isFaceInside(
      const int faceIndex,
      exahype::Vertex* const verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);

  /**
   * TODO(Dominic): Revise this docu.
   * Here we reset helper variables that play a role in
   * the neighbour merge methods.
   * These are the cell description attributes
   * riemannSolvePerformed[], and
   * faceDataExchangeCounter[].
   *
   * <h2>Shared Memory</h2>
   * The flag riemannSolvePerformed
   * indicates for every thread that touches a
   * face of a cell description if a Riemann Solve
   * was already performed for this face.
   *
   * <h2>MPI</h2>
   * This method resets Face data exchange counters:
   * To this end, we count the listings of a remote rank on each
   * of the faces surrounding a cell description.
   *
   * We perform the following actions depending on the counter value:
   * 0 - no connection: no send. Set to unreachable value.
   * 2^{d-2} - full face connection where cell is inside but half of face vertices are outside:
   * send at time of 2^{d-2}-th touch of face.
   * 4^{d-2} - full face connection where cell is inside and face vertices are all inside:
   * send at time of 2^{d-2}-th touch of face.
   * We require that #ifdef vertex.isBoundary() is set.
   */
  template <class CellDescription>
  static void resetNeighbourMergeHelperVariables(
      CellDescription& cellDescription,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
    for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
      cellDescription.setRiemannSolvePerformed(faceIndex,false);

      #ifdef Parallel
      int listingsOfRemoteRank =
          countListingsOfRemoteRankByInsideVerticesAtFace(
              faceIndex,fineGridVertices,fineGridVerticesEnumerator);
      if (listingsOfRemoteRank==0) {
        listingsOfRemoteRank = TWO_POWER_D;
      }
      cellDescription.setFaceDataExchangeCounter(faceIndex,listingsOfRemoteRank);
      assertion(cellDescription.getFaceDataExchangeCounter(faceIndex)>0);
      #endif
    }
  }

  template <class CellDescription>
  static void determineInsideAndOutsideFaces(
      CellDescription& cellDescription,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
    for (int faceIndex=0; faceIndex<DIMENSIONS_TIMES_TWO; faceIndex++) {
      // TODO(Dominic): Normally, isFaceInside has not to be called everytime here
      // but only once when the cell is initialised. Problem: Call addNewCellDescr.. from  merge..DueToForkOrJoin(...),
      // where no vertices are given. [Solved] - We send out the cellDescriptions from
      // the ADERDG/FV cell descr. heaps.
      cellDescription.setIsInside(faceIndex,isFaceInside(faceIndex,fineGridVertices,fineGridVerticesEnumerator));
    }
  }

  #ifdef Parallel
  /**
   * Returns true if the cell corresponding
   * to the vertices \p verticesAroundCell
   * is neighbour to a remote rank
   * via one of the faces.
   */
  static bool isAdjacentToRemoteRank(
      exahype::Vertex* const verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);

  /**
   * Count the listings of remote ranks sharing a vertex
   * adjacent to the face \p faceIndex of a cell with this rank.
   * This value is either 0, or 2^{d-1}.
   *
   * If we count 2^{d-1} listings, we directly know that this rank
   * shares a whole face with a remote rank.
   * If we count 0 listings, we do not have
   * a remote rank adjacent to this face.
   *
   * <h2>MPI</h2>
   * We know from the result of this function how
   * many vertices will try to exchange neighbour information
   * for this face.
   *
   * @developers:
   * TODO(Dominic): We currently check for uniqueness of the
   * remote rank. This might however not be necessary.
   */
  static int countListingsOfRemoteRankByInsideVerticesAtFace(
      const int faceIndex,
      exahype::Vertex* const verticesAroundCell,
      const peano::grid::VertexEnumerator& verticesEnumerator);
  #endif

  /**
   * Returns meta data describing the surrounding cell descriptions. The
   * routine is notably used by the automated adapters to derive adjacency
   * information on the cell level.
   */
  int getCellDescriptionsIndex() const;

  /**
   * TODO(Dominic): Add docu.
   */
  void setCellDescriptionsIndex(int cellDescriptionsIndex);

  /**
   * Returns the number of ADERDGCellDescriptions associated
   * with this cell.
   */
  int getNumberOfADERDGCellDescriptions() const;

  /**
   * Returns the number of FiniteVolumesCellDescriptions associated
   * with this cell.
   */
  int getNumberOfFiniteVolumeCellDescriptions() const;

  /**
   * Each cell points to a series of cell descriptions. The array holding the
   * series has to be stored on the heap and, consequently, initialised
   * properly. This is done by create() while destroy() cleans up. Please note
   * that you have to invoke create() once before you do anything with the cell
   * at all. You should destroy() in return in the very end.
   *
   * The operation shows that each cell in the tree can theoretically hold a
   * solver though only few do.
   *
   * This operation is used by addNewCellDescription().
   */
  void setupMetaData();

  /**
   * @see setupMetaData()
   */
  void shutdownMetaData();

  /**
   * todo docu
   *
   * setupMetaData() is called if cell hasn't been properly initialised before.
   */
  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::ADERDGCellDescription::Type cellType,
      const exahype::records::ADERDGCellDescription::RefinementEvent
          refinementEvent,
      const int level,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const tarch::la::Vector<DIMENSIONS, double>& cellOffset);

  void addNewCellDescription(
      const int solverNumber,
      const exahype::records::FiniteVolumesCellDescription::Type cellType,
//      const exahype::records::FiniteVolumesCellDescription::RefinementEvent
//          refinementEvent,    @todo Dominic
      const int level,
      const int parentIndex,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const tarch::la::Vector<DIMENSIONS, double>& cellOffset);

  /**
   * @return if this cell is initialised.
   *
   * @developers:
   * Note that it is not simply sufficient to
   * check if the heap index equals
   * multiscalelinkedcell::HangingVertexBookkeeper.
   *
   * @TODO bug
   */
  bool isInitialised() const;
};

#endif
