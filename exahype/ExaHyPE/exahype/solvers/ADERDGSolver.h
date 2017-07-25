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
 *
 * \author Dominic E. Charrier, Tobias Weinzierl, Fabian Güra
 **/

#ifndef _EXAHYPE_SOLVERS_ADERDG_SOLVER_H_
#define _EXAHYPE_SOLVERS_ADERDG_SOLVER_H_

#include <iostream>
#include <string>
#include <vector>

#include "exahype/solvers/Solver.h"

#include "peano/heap/Heap.h"
#include "peano/utils/Globals.h"

#include "tarch/Assertions.h"
#include "tarch/la/Vector.h"

#include "exahype/profilers/simple/NoOpProfiler.h"
#include "exahype/records/ADERDGCellDescription.h"

namespace exahype {
  namespace solvers {
    class ADERDGSolver;
  }
}

/**
 * Describes one solver.
 */
class exahype::solvers::ADERDGSolver : public exahype::solvers::Solver {
public:
  /**
   * Set to 0 if no floating point compression is used.
   */
  static double CompressionAccuracy;

  static bool SpawnCompressionAsBackgroundThread;

  #ifdef Asserts
  static double PipedUncompressedBytes;
  static double PipedCompressedBytes;
  #endif

  typedef exahype::DataHeap DataHeap;

  /**
   * Rank-local heap that stores ADERDGCellDescription instances.
   *
   * \note This heap might be shared by multiple ADERDGSolver instances
     that differ in their solver number and other attributes.
     @see solvers::Solver::RegisteredSolvers.
   */
  typedef exahype::records::ADERDGCellDescription CellDescription;
  typedef peano::heap::PlainHeap<CellDescription> Heap;

private:

  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * The number of unknowns/basis functions associated with each face of an
   * element.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerFace;

  /**
   * The total number of unknowns/basis functions associated with the 2^d faces
   * of an element.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerCellBoundary;

  /**
   * The total number of unknowns/basis functions associated with the volume of
   * a cell.
   * This number includes the unknowns of all state variables.
   */
  const int _unknownsPerCell;

  /**
   * The total number of volume flux unknowns/basis functions PLUS the number of
   * source unknowns
   * associated with the volume of a cell.
   * This number includes the unknowns of all state variables.
   *
   *
   */
  const int _fluxUnknownsPerCell;

  /**
   * The total number of space-time unknowns/basis functions associated with the
   * space-time volume of a cell and its time stepping interval.
   * This number includes the unknowns of all state variables.
   */
  const int _spaceTimeUnknownsPerCell;

  /**
   * The total number of space-time volume flux unknowns/basis functions
   * PLUS the number of space-time source unknowns associated with the
   * space-time volume of a cell and its time stepping interval.
   * This number includes the unknowns of all state variables.
   */
  const int _spaceTimeFluxUnknownsPerCell;

  /**
   * The size of data required to store cell volume based unknowns and
   * associated parameters.
   */
  const int _dataPerCell;

  /**
   * Minimum corrector time step size of all
   * cell descriptions in the previous iteration.
   */
  double _previousMinCorrectorTimeStepSize;

  /**
   * Minimum corrector time stamp of all cell descriptions.
   */
  double _minCorrectorTimeStamp;

  /**
   * Minimum predictor time stamp of all cell descriptions.
   * Always equal or larger than the minimum corrector time stamp.
   */
  double _minPredictorTimeStamp;

  /**
   * Minimum corrector time step size of
   * all cell descriptions.
   */
  double _minCorrectorTimeStepSize;

  /**
   * Minimum predictor time step size of
   * all cell descriptions.
   */
  double _minPredictorTimeStepSize;

  /**
   * Minimum next predictor time step size of
   * all cell descriptions.
   */
  double _minNextPredictorTimeStepSize;

  void tearApart(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa);
  void glueTogether(int numberOfEntries, int normalHeapIndex, int compressedHeapIndex, int bytesForMantissa);

  /**
   * Different to compress(), this operation is called automatically by
   * mergeNeighbours(). Therefore the routine is private.
   */
  void uncompress(exahype::records::ADERDGCellDescription& cellDescription);

  /**
   * TODO(Dominic): Add more docu.
   *
   * Mark a cell description of Cell for refinement or erasing based
   * on a user supplied physics based refinement criterion.
   *
   * <h2>Erasing</h2>
   * Note that we use a not so obvious strategy for performing
   * erasing operations. We first set an erasing request on
   * a parent cell description of type Ancestor or EmptyAncestor,
   * and then let its children of type Cell veto
   * this request if they want to keep their
   * solution or refine even further.
   *
   * We further veto erasing events if
   * a child of the parent itself is a parent
   * of cell descriptions of type Descendant/EmptyDescendant.
   *
   * <h2>Augmentation</h2>
   * Note that a cell description of type Cell is allowed to overwrite an augmentation request
   * by a refinement request if applicable.
   * The refinement event of a cell description of type Cell might be set to
   * an augmentation request in the methods mergeWithNeighbourData(...)
   * as well as in markForAugmentation(...) which is called from within
   * enterCell(...)
   */
  bool markForRefinement(
      CellDescription& pFine);

  /**
   * Check if cell descriptions of type Ancestor or Descendant need to hold
   * data or not based on virtual refinement criterion.
   * Then, allocate the necessary memory or deallocate the unnecessary memory.
   */
  void ensureOnlyNecessaryMemoryIsAllocated(
      CellDescription& fineGridCellDescription,
      const exahype::solvers::Solver::AugmentationControl& augmentationControl,
      const bool onMasterWorkerBoundary);

  /**
   * TODO(Dominic): Add docu.
   */
  bool markForAugmentation(
      CellDescription& pFine,
      const tarch::la::Vector<THREE_POWER_D, int>& neighbourCellDescriptionIndices,
      const bool onMasterWorkerBoundary);

  /*
   * Change the erasing children request to a change children to descendants
   * request of the coarse grid cell description's parent
   * if the coarse grid cell has children itself (of type Descendant).
   * Rationale: We cannot directly erase a Cell that has children (of type Descendant).
   *
   * Further, reset the deaugmenting children request if a coarse grid
   * Descendant has children (of type Descendant). Rationale:
   * We cannot erase a coarse grid cell that has children (of type Descendant)
   * before erasing the children.
   *
   * \note This operation spans over three spacetree levels. Calling
   * it requires that a cell description for
   * the same solver the \p coarseGridCellDescription is associated with
   * is registered on the fine grid cell.
   *
   * \note A more sophisticated procedure has to performed for the refinement event
   * AugmentationRequested. We need to use the taversal's descend event to handle
   * this event. We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
   * to check if we need to reset the deaugmenting request.
   *
   * TODO(Dominic): Make template function as soon as verified.
   */
  void vetoErasingOrDeaugmentingChildrenRequest(
      CellDescription& coarseGridCellDescription,
      const int fineGridCellDescriptionsIndex);

  /*
   * Resets the refinement event of a fine grid cell of type
   * Descendant to None if it was set to Refining.
   * The latter event indicates that the fine grid cells in
   * the next finer level have all been initialised with
   * type EmptyAncestor/Ancestor.
   *
   * Resets the augmentation event of a fine grid cell of type
   * Descendant to None if it was set to Augmenting.
   * The latter event indicates that the fine grid cells in
   * the next finer level have all been initialised with
   * type Descendant.
   *
   * TODO(Dominic): Erasing
   * TODO(Dominic): Make template function as soon as verified.
   */
  void startOrFinishCollectiveRefinementOperations(
      CellDescription& fineGridCellDescription);

  bool eraseCellDescriptionIfNecessary(
      const int cellDescriptionsIndex,
      const int fineGridCellElement,
      const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell,
      CellDescription& coarseGridCellDescription);

  /**
   * Initialise cell description of type Cell.
   * Initialise the refinement event with None.
   */
  void addNewCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const int coarseGridCellDescriptionsIndex,
      const int solverNumber);

  /**
   * Initialises helper cell descriptions of type Descendant
   * on the fine level after cell descriptions on the coarse level
   * have been flagged for augmentation and Peano has
   * created the requested new cells.
   *
   * Further sets the refinement event on a coarse grid Descendant to Augmenting
   * if the first new Descendant was initialised on the fine grid.
   */
  void addNewDescendantIfAugmentingRequested(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      CellDescription& coarseGridCellDescription,
      const int coarseGridCellDescriptionsIndex);

  /**
   * Initialises compute cell descriptions on the fine level (cell description type is Cell)
   * after coarse grid cell descriptions have been flagged for refinement and Peano has
   * created the requested new cells.
   * Erasing is not performed on cells belonging to the regular initial grid
   * of the solvers (see RegularMesh).
   *
   * Further sets the refinement event on a coarse grid Cell to Refining
   * if the first new Cell was initialised on the fine grid.
   */
  void addNewCellIfRefinementRequested(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const tarch::la::Vector<DIMENSIONS,int>& fineGridPositionOfCell,
      CellDescription& coarseGridCellDescription,
      const int coarseGridCellDescriptionsIndex);

  /**
   * Prolongates Volume data from a parent cell description to
   * \p cellDescription if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This method makes only sense for virtual shells
   * in the current AMR concept.
   */
  void prolongateVolumeData(
      CellDescription&       fineGridCellDescription,
      const CellDescription& coarseGridCellDescription,
      const tarch::la::Vector<DIMENSIONS, int>&      subcellIndex);

  /**
   * Restricts Volume data from \p cellDescriptio to
   * a parent cell description if the fine grid cell associated with
   * \p cellDescription is adjacent to a boundary of the
   * coarse grid cell associated with the parent cell description.
   *
   * \note This method makes only sense for real cells.
   * in the current AMR concept.
   */
  void restrictVolumeData(
      CellDescription&       coarseGridCellDescription,
      const CellDescription& fineGridCellDescription,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex);

  /**
   * Sets the face unknowns of a cell description of type Ancestor to zero.
   * This is typically done before we perform a face unknowns
   * restriction operation.
   */
  void prepareFaceDataOfAncestor(CellDescription& cellDescription);

  /**
   * Determine if the cell description of type
   * Descendant is on the cell boundary of its parent
   * of type Cell or Descendant with at least one of
   * its faces. If so restrict face data from the parent down
   * to the Descendant for those face(s).
   */
  void prolongateFaceDataToDescendant(
      CellDescription& cellDescription,
      SubcellPosition& SubcellPosition);

  /**
   * Solve the Riemann problem at the interface between two cells ("left" and
   * "right"). This method only performs a Riemann solve if at least one of the
   * cell descriptions (per solver) associated with the two cells is of type
   * ::Cell and none of the two cells belongs to the boundary.
   * In case a Riemann problem is solved,
   * the method further sets the ::riemannSolvePerformed
   * flags for the particular faces on both cell descriptions (per solver).
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * <h2>Rationale</h2>
   *
   * We did originally split up the boundary condition handling and the Riemann
   * updates into two mappings. This offers a functional decomposition. However,
   * both mappings then need a significiant number of technical administrative
   * code (cmp all the loops in touchVertexFirstTime and the redundant code to
   * manage the semaphores). We thus decided to merge both aspects. This also
   * should make sense from a performance point of view.
   *
   * We could potentially remove the face indices here if we had normals that
   * point outwards. However, we don't evaluate the direction of the normal and
   * thus need these counters as a Riemann problem on a face either could be
   * triggered by the left cell or by the right cell.
   *
   * \note The current implementation might classify cells with vertices that
   * are part of the
   * boundary of the domain or outside to be classified as inside of the domain
   * (volume-ratio based).
   *
   * \note We cannot solely check for indices of value
   * multiscalelinked::HangingVertexBookkepper::DomainBoundaryAdjacencyIndex
   * in vertex.getCellDescriptions() to determine if we are on the boundary of
   * the domain
   * since these values are overwritten by
   * multiscalelinked::HangingVertexBookkepper::RemoteAdjacencyIndex
   * if the domain boundary aligns with an MPI boundary
   * (see
   * multiscalelinkedcell::HangingVertexBookkeeper::updateCellIndicesInMergeWithNeighbour(...)).
   *
   * \note Not thread-safe.
   */
  void solveRiemannProblemAtInterface(
      CellDescription& pLeft,
      CellDescription& pRight,
      const int faceIndexLeft,
      const int faceIndexRight,
      double**  tempFaceUnknowns,
      double**  tempStateSizedVectors,
      double**  tempStateSizedSquareMatrices);

  /**
   * Apply the boundary conditions at the face with index \p faceIndex.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   *
   * @param[in] cellDescription         The cell description
   * @param[in] faceIndex               The index of the interface
   *                                    from the perspective of the cell/cell
   *                                    description. The index is computed as 2 times the
   *                                    position of the normal vector non-zero plus a
   *                                    value that encodes the normal vector direction
   *                                    (0 for negative direction, 1 for positive direction).
   * \note Not thread-safe.
   */
  void applyBoundaryConditions(
      CellDescription& p,
      const int faceIndex,
      double**  tempFaceUnknowns,
      double**  tempStateSizedVectors,
      double**  tempStateSizedSquareMatrices);

  /**
   * Checks if no unnecessary memory is allocated for the cell description.
   * If this is not the case, it deallocates the unnecessarily allocated memory.
   *
   * This operation is thread safe as we serialise it.
   */
  void ensureNoUnnecessaryMemoryIsAllocated(CellDescription& cellDescription);

  /**
   * Checks if all the necessary memory is allocated for the cell description.
   * If this is not the case, it allocates the necessary
   * memory for the cell description.
   */
  void ensureNecessaryMemoryIsAllocated(exahype::records::ADERDGCellDescription& cellDescription);

#ifdef Parallel
  /**
   * Data messages per neighbour communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerNeighbourCommunication;
  /**
   * Data messages per fork/join communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerForkOrJoinCommunication;
  /**
   * Data messages per master worker communication.
   * This information is required by the sendEmpty...(...)
   * methods.
   */
  static const int DataMessagesPerMasterWorkerCommunication;

  /**
   * Single-sided version of the other solveRiemannProblemAtInterface(). It
   * works only on one cell and one solver within this cell and in return
   * hands in the F and Q values explicitly through  indexOfQValues and
   * indexOfFValues. The Riemann solver is invoked and the bits are set
   * accordingly no matter of what they did hold before, i.e. different to
   * the standard solveRiemannProblemAtInterface() operation, we do not
   * check whether we shall run a Riemann solver or not.
   *
   * This method further synchronises the ADERDGCellDescription
   * with the corresponding solver if this is required by the time stepping
   * scheme.
   * This operation must be performed in mergeWithNeighbour(...) and
   * touchVertexFirstTime(...) since both callbacks touch the
   * ADERDGCellDescriptions before the other callbacks.
   *
   * \note Not thread-safe.
   */
  void solveRiemannProblemAtInterface(
      records::ADERDGCellDescription& cellDescription,
      const int faceIndex,
      const int indexOfQValues,
      const int indexOfFValues,
      double**  tempFaceUnknowns,
      double**  tempStateSizedVectors,
      double**  tempStateSizedSquareMatrices);

  /**
   * Sets heap indices of all ADER-DG cell descriptions that were
   * received due to a fork or join event to
   * multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
   * and the parent index of the cell descriptions to the specified \p
   * parentIndex.
   */
  static void resetDataHeapIndices(
      const int cellDescriptionsIndex,
      const int parentIndex);

  /**
   * Checks if the parent index of a fine grid cell description
   * was set to RemoteAdjacencyIndex during a previous forking event.
   *
   * If so, check if there exists a coarse grid cell description
   * which must have been also received during a previous fork event.
   * If so, update the parent index of the fine grid cell description
   * with the coarse grid cell descriptions index.
   */
  void ensureConsistencyOfParentIndex(
      CellDescription& cellDescription,
      const int coarseGridCellDescriptionsIndex,
      const int solverNumber);

#endif
  /**
   * Run over the persistent fields of the ADER-DG cell and determine the
   * average per unknown.' The result is stored within
   *
   */
  void determineUnknownAverages(exahype::records::ADERDGCellDescription& cellDescription);

  /**
   * Runs over all entries and adds sign times the average value. So if you
   * hand in a -1, you compute the hierarchical transform. If you hand in a +1,
   * you compute the inverse hierarchical transform.
   */
  void computeHierarchicalTransform(exahype::records::ADERDGCellDescription& cellDescription, double sign);

  /**
   * This routine runs over the unknowns, asks the Heap's compression routines
   * to identify a reasonable compression level, stores this value and
   * afterwrads pipes the dofs into the byte stream. If you don't run with
   * assertions, the code does clear the double heap data afterwards. With
   * assertions, we leave it there and thus allow pullUnknownsFromByteStream()
   * to do quite some validation.
   */
  void putUnknownsIntoByteStream(exahype::records::ADERDGCellDescription& cellDescription);

  /**
   *
   *
   * <h2>Multicore</h2>
   *
   * Unknowns are pulled from the input stream indirectly through
   * touchVertexFirstTime(). It always recycles heap data, so can be triggered
   * while other routines already do something with the cell description. There
   * usually are enough entries to recycle available.
   *
   * However, it may happen that we run out of recycled entries if we run into
   * a large regular subgrid for the first time. We can identify a run out as
   * we get a -1 from the heap. In this case, there are a couple of things to
   * do.
   *
   * - Wait for any background task to finish. Other parts of the grid might
   *   have triggered a compression in the background. So we have to wait for
   *   those guys to finish, as they rely on an invariant heap.
   * - Lock very pessimistically. No two operations (only touchVertexFirstTime
   *   calls should run in parallel, but I'm not 100% sure) should run.
   * - Create additional data.
   */
  void pullUnknownsFromByteStream(exahype::records::ADERDGCellDescription& cellDescription);

  class CompressionTask {
    private:
      ADERDGSolver&                             _solver;
      exahype::records::ADERDGCellDescription&  _cellDescription;
    public:
      CompressionTask(
        ADERDGSolver&                             _solver,
        exahype::records::ADERDGCellDescription&  _cellDescription
      );

      void operator()();
  };

public:
  /**
   * Returns the ADERDGCellDescription.
   */
  static Heap::HeapEntries& getCellDescriptions(
      const int cellDescriptionsIndex) {
    assertion1(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex);

    return Heap::getInstance().getData(cellDescriptionsIndex);
  }

  /**
   * Returns the ADERDGCellDescription.
   */
  static CellDescription& getCellDescription(
      const int cellDescriptionsIndex,
      const int element) {
    assertion2(Heap::getInstance().isValidIndex(cellDescriptionsIndex),cellDescriptionsIndex,element);
    assertion2(element>=0,cellDescriptionsIndex,element);
    assertion2(static_cast<unsigned int>(element)<Heap::getInstance().getData(cellDescriptionsIndex).size(),cellDescriptionsIndex,element);

    return Heap::getInstance().getData(cellDescriptionsIndex)[element];
  }

  /**
   * Returns if a ADERDGCellDescription type holds face data.
   */
  static bool holdsFaceData(const CellDescription::Type& cellDescriptionType) {
    return cellDescriptionType==CellDescription::Cell       ||
        cellDescriptionType==CellDescription::Ancestor   ||
        cellDescriptionType==CellDescription::Descendant;
  }

  ADERDGSolver(
      const std::string& identifier,
      int numberOfVariables, int numberOfParameters, int nodesPerCoordinateAxis,
      double maximumMeshSize,
      exahype::solvers::Solver::TimeStepping timeStepping,
      std::unique_ptr<profilers::Profiler> profiler =
          std::unique_ptr<profilers::Profiler>(
              new profilers::simple::NoOpProfiler("")));

  virtual ~ADERDGSolver() {}

  // Disallow copy and assignment
  ADERDGSolver(const ADERDGSolver& other) = delete;
  ADERDGSolver& operator=(const ADERDGSolver& other) = delete;
  /**
   * Returns the maximum extent a mesh cell is allowed to have
   * in all coordinate directions.
   * This maximum mesh size is used both as a
   * constraint on the AMR as well as to set up the initial
   * grid. If you return the extent of the computational domain in
   * each coordinate direction or larger values,
   * you indicate that this solver is not active in the domain.
   */
  double getMaximumMeshSize() const;

  /**
   * Returns the number of state variables.
   */
  int getNumberOfVariables() const;

  /**
   * Returns the number of parameters, e.g.,material constants etc.
   */
  int getNumberOfParameters() const;

  /**
   * This operation returns the number of space time
   * unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeUnknownsPerCell() const;

  /**
   * This operation returns the number of space time
   * flux unknowns per cell.
   *
   * Note that this operation might only have a meaning for space-time type
   * discretisation methods.
   */
  int getSpaceTimeFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns per cell located in
   * the interior of a cell.
   */
  int getUnknownsPerCell() const;

  /**
   * This operation returns the number of flux unknowns per cell
   * located in the interior of a cell.
   */
  int getFluxUnknownsPerCell() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of the boundary of a cell.
   */
  int getUnknownsPerCellBoundary() const;

  /**
   * This operation returns the number of unknowns that are located
   * on or in the vicinity of each face of a cell.
   */
  int getUnknownsPerFace() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  /**
   * This operation returns the size of data required
   * to store cell volume based unknowns and associated parameters.
   */
  int getDataPerCell() const;
  
  /**
   * Getter for the size of the array allocated that can be overriden
   * to change the allocated size independently of the solver parameters.
   * For example to add padding forthe optimised kernel
   */
  virtual int getTempSpaceTimeUnknownsSize()     const {return getSpaceTimeUnknownsPerCell()+getUnknownsPerCell();}
  virtual int getTempSpaceTimeFluxUnknownsSize() const {return getSpaceTimeFluxUnknownsPerCell();}
  virtual int getTempUnknownsSize()              const {return getUnknownsPerCell();}
  virtual int getTempFluxUnknownsSize()          const {return getFluxUnknownsPerCell();}
  virtual int getBndFaceSize()                   const {return getUnknownsPerFace();}
  virtual int getBndTotalSize()                  const {return getUnknownsPerCellBoundary();}
  
  virtual bool alignTempArray()                  const {return false;}

  /**
   * False for generic solver, may be true for optimized one
   * Used only for debug assertions
   */
  virtual bool usePaddedData_nVar() const {return false;}
  virtual bool usePaddedData_nDoF() const {return false;}
  
  //TODO KD
  virtual bool isDummyKRequired() const {return false;}
  
  /**
   * @brief Adds the solution update to the solution.
   *
   * @param[inout] luh  Cell-local solution DoF.
   * @param[in]    lduh Cell-local update DoF.
   * @param[dt]    dt   Time step size.
   */
  virtual void solutionUpdate(double* luh, const double* const lduh,
                              const double dt) = 0;

  /**
   * @brief Computes the volume flux contribution to the cell update.
   *
   * @param[inout] lduh Cell-local update DoF.
   * @param[in]    cellSize   Extent of the cell in each coordinate direction.
   * @param[dt]    dt   Time step size.
   */
  virtual void volumeIntegral(
      double* lduh, const double* const lFhi,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) = 0;

  /**
   * @brief Computes the surface integral contributions
   * to the cell update.
   *
   * @param[inout] lduh   Cell-local update DoF.
   * @param[in]    lFhbnd Cell-local DoF of the boundary extrapolated fluxes.
   * @param[in]    cellSize     Extent of the cell in each coordinate direction.
   */
  virtual void surfaceIntegral(
      double* lduh, const double* const lFhbnd,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) = 0;

  /**
   * @brief Computes the normal fluxes (or fluctuations) at the interface of two
   *cells.
   *
   * @param[inout] FL             Flux DoF belonging to the left cell.
   * @param[inout] FR             Flux DoF belonging the right cell.
   * @param[in]    QL             DoF of the boundary extrapolated predictor
   *                              belonging to the left cell.
   * @param[in]    QR             DoF of the boundary extrapolated predictor
   *                              belonging to the right cell.
   * @param[in]    tempFaceUnknownsArray        Temporary array of the size of a face unknowns array.
   * @param[in]    tempStateSizedVectors        Five (5) state sized (=number of variables) temporary variables.
   * @param[in]    tempStateSizedSquareMatrices Three (3) temporary variables of the size number of variables squared.
   * @param[in]    normalNonZero  Index of the nonzero normal vector component,
   *i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   */
  virtual void riemannSolver(double* FL, double* FR, const double* const QL,
                             const double* const QR,
                             double*   tempFaceUnknownsArray,
                             double**  tempStateSizedVectors,
                             double**  tempStateSizedSquareMatrices,
                             const double dt,
                             const int normalNonZero) = 0;

  /**
   * Return the normal fluxes (or fluctuations) and state variables at the boundary.
   *
   * @param[inout] fluxOut       Flux DoF belonging to the left cell.
   * @param[inout] stateOut      DoF of the boundary extrapolated predictor
   *                             belonging to the left cell.
     @param[in]    fluxIn        Flux DoF belonging to the left cell.
   * @param[in]    stateIn       DoF of the boundary extrapolated predictor
   *                             belonging to the left cell.
   * @param[in]    cellCentre    Cell centre.
   * @param[in]    cellSize      Cell size.
   * @param[in]    t             The time.
   * @param[in]    dt            A time step size.
   * @param[in]    normalNonZero Index of the nonzero normal vector component,
   *i.e., 0 for e_x, 1 for e_y, and 2 for e_z.
   */
  virtual void boundaryConditions(double* fluxOut,
                                  double* stateOut,
                                  const double* const fluxIn,
                                  const double* const stateIn,
                                  const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
                                  const tarch::la::Vector<DIMENSIONS,
                                  double>& cellSize,
                                  const double t,const double dt,
                                  const int faceIndex,
                                  const int normalNonZero) = 0;


  /**
   * @brief Computes cell-local predictor space-time, volume, and face DoF.
   *
   * Computes the cell-local space-time predictor lQi, the space-time volume
   *flux lFi,
   * the predictor lQhi, the volume flux lFhi, the boundary
   * extrapolated predictor lQhbnd and normal flux lFhbnd.
   *
   * @param[inout] lQi       Space-time predictor DoF.
   * @param[in]    lQi_old   Old space-time predictor DoF - only used in Picard loop.
   * @param[in]    rhs       The right-hand side vector - only used in Picard loop.
   * @param[in]    rhs_0     Constant term of the right-hand side - only used in Picard loop.
   * @param[inout] lFi       Space-time flux DoF.
   * @param[inout] lQhi      Predictor DoF
   * @param[inout] lFhi      Volume flux DoF.
   * @param[out]   luh       Solution DoF.
   * @param[in]    cellSize     Extent of the cell in each coordinate direction.
   * @param[in]    dt     Time step size.
   */
  virtual void spaceTimePredictor(
      double*  lQhbnd, double* lFhbnd,
      double** tempSpaceTimeUnknowns,
      double** tempSpaceTimeFluxUnknowns,
      double*  tempUnknowns,
      double*  tempFluxUnknowns,
      double*  tempStateSizedVector,
      const double* const luh,
      const tarch::la::Vector<DIMENSIONS, 
      double>& cellSize, 
      const double dt,
      double* pointForceSources) = 0;

  /**
   * \brief Returns a stable time step size.
   *
   * \param[in] luh             Cell-local solution DoF.
   * \param[in] tempEigenvalues A temporary array of size equalling the number of variables.
   * \param[in] cellSize        Extent of the cell in each coordinate direction.
   */
  virtual double stableTimeStepSize(
      const double* const luh,
      double* tempEigenvalues,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize) = 0;

  /**
   * This operation allows you to impose time-dependent solution values
   * as well as to add contributions of source terms.
   * Please be aware that this operation is called per time step if
   * the corresponding predicate hasToUpdateSolution() yields true for the
   * region and time interval.
   *
   * \param t  The new time stamp after the solution update.
   * \param dt The time step size that was used to update the solution.
   *           This time step size was computed based on the old solution.
   *           If we impose initial conditions, i.e, t=0, this value
   *           equals std::numeric_limits<double>::max().
   */
  virtual void solutionAdjustment(
      double* luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt) = 0;

  /**
   * This hook can be used to trigger solution adjustments within the
   * region corresponding to \p cellCentre and \p dx
   * and the time interval corresponding to t and dt.
   *
   * \param t  The new time stamp after the solution update.
   * \param dt The time step size that was used to update the solution.
   *           This time step size was computed based on the old solution.
   *           If we impose initial conditions, i.e, t=0, this value
   *           equals std::numeric_limits<double>::max().
   */
  virtual bool hasToAdjustSolution(
      const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& dx,
      const double t,
      const double dt) = 0;

  /**
   * @defgroup AMR Solver routines for adaptive mesh refinement
   */
  ///@{
  /**
   * The refinement criterion that must be defined by the user.
   *
   */
  // @todo: 16/04/06:Dominic Etienne Charrier Consider to correct the level in
  // the invoking code, i.e., level-> level-1
  // since this is was the user expects.
  virtual exahype::solvers::Solver::RefinementControl refinementCriterion(
      const double* luh, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
      const tarch::la::Vector<DIMENSIONS, double>& cellSize,
      const double time,
      const int level) = 0;

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels can
   * be larger than one. Let \f$l\f$ be the level difference. The
   * vector \p subfaceIndex does contain values in the range
   * \f$0,1,\ldots,3^l-1\f$.
   */
  virtual void faceUnknownsProlongation(
      double* lQhbndFine, double* lFhbndFine, const double* lQhbndCoarse,
      const double* lFhbndCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex) = 0;

  /**
   * Restricts fine grid face unknowns on level \p fineGridLevel
   * up to level \p coarseGridLevel and adds them to the coarse grid unknowns.
   *
   * \note For the considered AMR concept, the difference in levels can
   * be larger than one. Let \f$l\f$ be the level difference. The
   * vector \p subfaceIndex does contain values in the range
   * \f$0,1,\ldots,3^l-1\f$.
   */
  virtual void faceUnknownsRestriction(
      double* lQhbndCoarse, double* lFhbndCoarse, const double* lQhbndFine,
      const double* lFhbndFine, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS - 1, int>& subfaceIndex) = 0;

  /**
   * Project coarse grid face unknowns
   * on level \p coarseGridLevel down to level \p fineGridLevel
   * and writes them to the fine grid unknowns
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsProlongation(
      double* luhFine, const double* luhCoarse, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;

  /**
   * Restricts fine grid volume unknowns on level \p fineGridLevel
   * up to level \p coarseGridLevel and adds them to the coarse grid unknowns.
   *
   * \note For the considered AMR concept, the difference in levels is always
   * equal to one. The vector \p subcellIndex does contain values in the range
   * \f$0,1,2\f$.
   */
  virtual void volumeUnknownsRestriction(
      double* luhCoarse, const double* luhFine, const int coarseGridLevel,
      const int fineGridLevel,
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex) = 0;
  ///@}
  
  //TODO KD
  virtual void dummyK_GeneratedCall(
    const double t,
    const double dt, 
    const tarch::la::Vector<DIMENSIONS,double>& center,
    const tarch::la::Vector<DIMENSIONS,double>& dx, 
    double* tempPointForceSources);


  /**
   * A criterion determining if the degrees of freedoms of
   * the cell-wise solution luh are physically admissible.
   *
   * \note We require that the cell-local minimum and maximum
   * of the solution values has been computed
   * a-priori.
   *
   * This operation is required for limiting.
   */
  virtual bool physicalAdmissibilityDetection(const double* const QMin,const double* const QMax) { return true; }

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * meta data.
   */
  void synchroniseTimeStepping(
      CellDescription& p) const;

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * meta data.
   *
   * \param[in] element Index of the cell description in
   *                    the array at address \p cellDescriptionsIndex.
   */
  void synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) override;

  /**
   * Update and reset corrector and predictor
   * time stamps and time step sizes according to the chosen
   * time stepping variant.
   *
   * Further reset the minimum and maximum cell sizes
   * to numeric limit values.
   *
   * \note The minimum and maximum cell sizes do
   * not need to be reset to numeric limit values
   * in every time step for uniform mesh refinement
   * static adaptive mesh refinement
   * but we still do it since we want to
   * utilise dynamic adaptive mesh refinement
   * since we want to dynamic adaptive mesh refinement
   * eventually.
   */
  void startNewTimeStep() override;

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep();

  /**
   * Similar to reconstructStandardTimeSteppingData for roll backs.
   */
  void reconstructStandardTimeSteppingDataAfterRollback();

  /**
   * If we use the original time stepping
   * scheme with multiple algorithmic phases,
   * we have to call this method
   * after the time step computation phase.
   *
   * It ensures that the corrector time stamp
   * and step size equal their predictor
   * equivalents.
   */
  void reconstructStandardTimeSteppingData();

  /**
   * After the mesh has been updated,
   * reset the predictor time stamp to
   * the value corrector time stamp plus
   * nextPredictorTimeStepSize which is
   * the newly computed time step size.
   * Furthermore, set the predictor time
   * step size to nextPredictorTimeStepSize.
   *
   * This whole procedure is a little confusing but
   * necessary due to the shift
   * of the ADER-DG phases in our implementation.
   */
  void reinitialiseTimeStepData() override;

  /**
   * Update predictor time step size
   *
   * This operation takes the minimum of the current predictor time step size
   * and the argument handed in. The routine is used in
   * TimeStepComputation to determine the subsequent time step size.
   *
   * <h1>Thread-safety</h1>
   *
   * This operation is not thread safe.
   *
   */
  void updateMinNextPredictorTimeStepSize(
      const double& nextPredictorTimeStepSize);

  double getMinNextPredictorTimeStepSize() const;

  // todo 16/02/25:Dominic Etienne Charrier: It follows stuff that must be
  // revised:

  // todo 25/02/16:Dominic Etienne Charrier
  // Remove the time stamps that are not used in ExaHype.
  void setMinCorrectorTimeStamp(double minCorectorTimeStamp);

  double getMinCorrectorTimeStamp() const;

  void setMinPredictorTimeStamp(double minPredictorTimeStamp);

  double getMinPredictorTimeStamp() const;

  void setMinCorrectorTimeStepSize(double minCorrectorTimeStepSize);

  double getMinCorrectorTimeStepSize() const;

  void setMinPredictorTimeStepSize(double minPredictorTimeStepSize);

  double getMinPredictorTimeStepSize() const;

  double getPreviousMinCorrectorTimeStepSize() const;

  void setPreviousMinCorrectorTimeStepSize(double value);

  double getMinTimeStamp() const override {
    return getMinCorrectorTimeStamp();
  }

  double getMinTimeStepSize() const override {
    return getMinCorrectorTimeStepSize();
  }

  double getMinNextTimeStepSize() const override {
    return getMinNextPredictorTimeStepSize();
  }

  void updateMinNextTimeStepSize( double value ) override {
    updateMinNextPredictorTimeStepSize(value);
  }

  void initSolverTimeStepData(double value) override {
    setPreviousMinCorrectorTimeStepSize(0.0);
    setMinCorrectorTimeStepSize(0.0);
    setMinPredictorTimeStepSize(0.0);
    setMinCorrectorTimeStamp(value);
    setMinPredictorTimeStamp(value);
  }

  bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const override {
    return Heap::getInstance().isValidIndex(cellDescriptionsIndex);
  }

  /**
   * @todo Dominic, kannst Du mir reinschreiben, was die Routine tut und warum die so heisst? Der Name suggeriert, dass was schiefgehen kann
   */
  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override;

  SubcellPosition computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) override;

  ///////////////////////////////////
  // MODIFY CELL DESCRIPTION
  ///////////////////////////////////
  bool enterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  bool leaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) override;

  ///////////////////////////////////
  // CELL-LOCAL
  ///////////////////////////////////
  /**
   * Computes the space-time predictor quantities, extrapolates fluxes
   * and (space-time) predictor values to the boundary and
   * computes the volume integral.
   *
   * \param[in] tempSpaceTimeUnknows      Array of size 4 containing space-time predictor sized temporary arrays (see nonlinear predictor kernel).
   * \param[in] tempSpaceTimeFluxUnknowns Array of size 2 containing space-time predictor volume flux sized temporary arrays (see linear predictor kernel).
   * \param[in] tempUnknowns              Solution sized temporary array.
   * \param[in] tempFluxUnknowns          Volume flux sized temporary array.
   * \param[in] tempStateSizedVector      A vector of size of the state vector (=number of variables).
   */
  void performPredictionAndVolumeIntegral(
      exahype::records::ADERDGCellDescription& cellDescription,
      double** tempSpaceTimeUnknowns,
      double** tempSpaceTimeFluxUnknowns,
      double*  tempUnknowns,
      double*  tempFluxUnknowns,
      double*  tempStateSizedVector,
      double*  tempPointForceSources);

  void validateNoNansInADERDGSolver(
      const CellDescription& cellDescription,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      const std::string&                   methodTraceOfCaller);

  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) override;


  /**
   * If we use the original time stepping
   * scheme with multiple algorithmic phases,
   * we have to call this method
   * after a time step update with startNewTimeStep().
   *
   * It ensures that the corrector time size
   * is set to the admissible time step size that
   * was computed using the latest corrector solution.
   */
  void reconstructStandardTimeSteppingData(const int cellDescriptionsIndex,int element) const;

  /**
   * Rolls the solver time step data back to the
   * previous time step for a cell description.
   * Note that the newest time step
   * data is lost in this process.
   */
  void rollbackToPreviousTimeStep(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Similar to reconstructStandardTimeSteppingData for roll backs
   */
  void reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int element) const;

  /**
   * <h2>Solution adjustments</h2>
   * The solution is initially at time
   * cellDescription.getCorrectorTimeStamp().
   * The value cellDescription.getCorrectorTimeStepSize()
   * has initially no meaning and
   * equals std::numeric_limits<double>::max().
   */
  void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /*
   * Simply adds the update degrees of freedom
   * to the solution degrees of freedom.
   * Does not compute the surface integral.
   *
   * \deprecated We will not store the update field anymore
   * but a previous solution.
   */
  void addUpdateToSolution(
      CellDescription& cellDescription,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  /**
   * Computes the surface integral contributions to the
   * cell update and then adds the update degrees
   * on the solution degrees of freedom.
   *
   * <h2>Solution adjustments</h2>
   * After the update, the solution is at time
   * cellDescription.getCorrectorTimeStamp() + cellDescription.getCorrectorTimeStepSize().
   * The value cellDescription.getCorrectorTimeStepSize()
   * handed to the solution adjustment function is the one
   * used to update the solution.
   *
   * \todo We will not store the update field anymore
   * but a previous solution. We will thus only perform
   * a solution adjustment and adding of source term contributions here.
   */
  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * Rolls back the solver's solution on the
   * particular cell description.
   * This method is used by the ADER-DG a-posteriori
   * subcell limiter.
   *
   * Uses the corrector time step size to perform the rollback.
   * Thus make sure to invoke ::rollbackToPreviousTimeStepSize() beforehand
   * if the patch has already advanced to next time step.
   *
   * <h2>Open issues</h2>
   * A rollback is of course not possible if we have adjusted the solution
   * values. Assuming the rollback is invoked by a LimitingADERDGSolver,
   * we should use the adjusted FVM solution as reference solution.
   * A similar issue occurs if we impose initial conditions that
   * include a discontinuity.
   *
   * TODO(Dominic): A rollback is of course not possible if we have adjusted the solution
   * values. In this case, we should use the adjusted FVM solution as reference.
   * A similar issue occurs if we impose the initial conditions.
   */
  void rollbackSolution(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

  void swapSolutionAndPreviousSolution(
      const int cellDescriptionsIndex,
      const int element) const;

  void preProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void postProcess(
      const int cellDescriptionsIndex,
      const int element) override;

  void prolongateDataAndPrepareDataRestriction(
      const int cellDescriptionsIndex,
      const int element) override;

  void restrictData(
      const int cellDescriptionsIndex,
      const int element,
      const int parentCellDescriptionsIndex,
      const int parentElement,
      const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) override;

  ///////////////////////////////////
  // PARENT<->CHILD
  ///////////////////////////////////
  // TODO(Dominic): Extract prolongation and restriction operations.

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeNeighbours(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  void mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  void prepareNextNeighbourMerging(
      const int cellDescriptionsIndex,const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const override;
#ifdef Parallel
  /**
   * Sends all the cell descriptions at address \p
   * cellDescriptionsIndex to the rank \p toRank.
   *
   * <h2>Adaptive mesh refinement</h2>
   * For adaptive meshes, we further fix the type
   * of a descendant to RemoteBoundaryDescendant
   * at both sides of master-worker boundaries.
   *
   * We further fix the type of an Ancestor
   * to RemoteBoundaryAncestor if the parent
   * of the cell description on the master side
   * is also of type RemoteBoundaryAncestor or an
   * Ancestor.
   *
   * \note The data heap indices of the cell descriptions are not
   * valid anymore on rank \p toRank.
   */
  static void sendCellDescriptions(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Sends an empty message to the rank \p toRank.
   */
  static void sendEmptyCellDescriptions(
      const int                                    toRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Receives cell descriptions from rank \p fromRank
   * and resets the data heap indices to -1.
   *
   * If a received cell description has the same
   * solver number as a cell description in the
   * array at address \p cellDescriptionsIndex,
   * we merge the metadata (time stamps, time step size)
   * of both cell descriptions.
   *
   * If no cell description in the array at address
   * \p cellDescriptionsIndex can be found with the
   * same solver number than a received cell description,
   * we push the received cell description to
   * the back of the array at address \p cellDescriptions
   * Index.
   *
   * This operation is intended to be used in combination
   * with the solver method mergeWithWorkerOrMasterDataDueToForkOrJoin(...).
   * Here, we would merge first the cell descriptions sent by the master and worker
   * and then merge the data that is sent out right after.
   */
  static void mergeCellDescriptionsWithRemoteData(
      const int                                    fromRank,
      exahype::Cell&                               localCell,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  /**
   * Drop cell descriptions received from \p fromRank.
   */
  static void dropCellDescriptions(
      const int                                    fromRank,
      const peano::heap::MessageType&              messageType,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level);

  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeWithNeighbourMetadata(
      const int neighbourTypeAsInt,
      const int cellDescriptionsIndex,
      const int element) override;

  void sendDataToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendEmptyDataToNeighbour(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithNeighbourData(
      const int                                    fromRank,
      const int                                    neighbourTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      double**                                     tempFaceUnknowns,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  /**
   * Drop solver data from neighbour rank.
   */
  void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  ///////////////////////////////////
  // FORK OR JOIN
  ///////////////////////////////////

  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;


  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////
  bool hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) override;

  void sendDataToMaster(
      const int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerData(
      const int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToMaster(
      const int                                    masterRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendEmptyDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerData(
      const int                                     workerRank,
      const int                                     workerTypeAsInt,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////

  void sendDataToWorker(
      const                                        int workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithMasterData(
      const                                        int masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToWorker(
      const int                                    workerRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendEmptyDataToWorker(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithMasterData(
      const int                                    masterRank,
      const int                                    masterTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void dropMasterData(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;
#endif

  std::string toString() const override;

  void toString (std::ostream& out) const override;

  /**
   * The counterpart of uncompress.
   *
   * <h2> Shared memory parallelisation </h2>
   *
   * Different to the compression, we don't have to take care about any races:
   * the compression is invoked by enterCell or leaveCell respectively, i.e.
   * exactly once per cell. This can happen in parallel for multiple cells
   * which is fine.
   *
   * However, we have to take care about the interplay of compression and
   * uncompression. The
   */
  void compress(exahype::records::ADERDGCellDescription& cellDescription);
};

#endif
