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
 
#ifndef _EXAHYPE_SOLVERS_SOLVER_H_
#define _EXAHYPE_SOLVERS_SOLVER_H_

#include <memory>
#include <string>
#include <iostream>
#include <vector>

#include "tarch/la/Vector.h"
#include "tarch/multicore/BooleanSemaphore.h"

#include "peano/utils/Globals.h"
#include "peano/grid/VertexEnumerator.h"
#include "peano/heap/DoubleHeap.h"

#include "exahype/profilers/Profiler.h"
#include "exahype/profilers/simple/NoOpProfiler.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

// todo 08/02/16:Dominic Etienne Charrier
// move somewhere else
// is evaluated at compile time
constexpr int power(int basis, int exp) {
  return (exp == 0) ? 1 : basis * power(basis, exp - 1);
}

#ifdef ALIGNMENT
constexpr int addPadding(const int originalSize) {
  return ALIGNMENT/8 * static_cast<int>((originalSize+(ALIGNMENT/8-1))/(ALIGNMENT/8));
}
#else
constexpr int addPadding(const int originalSize) {
  return originalSize;
}
#endif

namespace exahype {
  // Forward declarations
  class Cell;
  class Vertex;
  /**
   * We store the degrees of freedom associated with the ADERDGCellDescription and FiniteVolumesCellDescription
   * instances on this heap.
   * We further use this heap to send and receive face data from one MPI rank to the other.
   */
  typedef peano::heap::PlainDoubleHeap DataHeap;

  namespace solvers {
    class Solver;

    typedef std::vector<Solver*> RegisteredSolversEntries;
    /**
     * All the registered solvers. Has to be declared extern in C++ standard as
     * it is instantiated in the corresponding cpp file.
     */
    // TODO: std::vector<std::unique_ptr<Solver>> ?!
    extern std::vector<Solver*> RegisteredSolvers;
}  // namespace solvers
}  // namespace exahype

/**
 * Describes one solver.
 */
class exahype::solvers::Solver {
 public:
  /**
   * Some solvers can deploy data conversion into the background. How this is
   * done is solver-specific. However, we have to wait until all tasks have
   * terminated if we want to modify the heap, i.e. insert new data or remove
   * data. Therefore, the wait (as well as the underlying semaphore) belong
   * into this abstract superclass.
   */
  static void waitUntilAllBackgroundTasksHaveTerminated();

  /**
   * The type of a solver.
   */
  enum class Type { ADERDG, FiniteVolumes, LimitingADERDG }; // TODO(Dominic): Get rid of the underscore

  /**
   * The time stepping mode.
   */
  enum class TimeStepping {
    /**
     * In the global time stepping mode, every cells works with the same time step.
     */
    Global,
    /**
     * In the fixed time stepping mode, we assume that each cell advanced in
     * time with the prescribed time step size. No CFL condition is checked.
     */
    GlobalFixed
    // Local, Anarchic
  };

  /**
   * The refinement control.
   */
  enum class RefinementControl { Keep = 0, Refine = 1, Erase = 2 };

  /**
   * This struct is used in the AMR context
   * to lookup a parent cell description and
   * for computing the subcell position of the child
   * with respect to this parent.
   *
   * TODO(Dominic): Move to more appropriate place.
   */
  typedef struct SubcellPosition {
    int parentCellDescriptionsIndex;
    int parentElement;
    tarch::la::Vector<DIMENSIONS,int> subcellIndex;
    int levelDifference;

    SubcellPosition() :
      parentCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex),
      parentElement(NotFound),
      subcellIndex(-1),
      levelDifference(-1) {}
    ~SubcellPosition() {}
  } SubcellPosition;

  /**
   * The augmentation control states.
   */
  enum class AugmentationControl {
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCell = 0,
    /**
     * Indicates that a spacetree cell is next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor.
     */
    NextToAncestor = 1,
    /**
     * Indicates that a spacetree cell is both, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, and
     * a spacetree cell of type
     * exahype::records::ADERDGCellDescription::Cell.
     */
    NextToCellAndAncestor = 2,
    /**
     * Indicates that a spacetree cell is neither, next to another spacetree cell
     * of type exahype::records::ADERDGCellDescription::Ancestor or
     * exahype::records::ADERDGCellDescription::EmptyAncestor, nor
     * next to a spacetree cell of type exahype::records::ADERDGCellDescription::Cell.
     *
     * A cell of type exahype::records::ADERDGCellDescription::Descendant can then request erasing.
     * A cell of type exahype::records::ADERDGCellDescription::Cell does then not need
     * to request augmenting.
     */
        Default = 3
  };

  /**
   * Default return value of function getElement(...)
   * If we do not find the element in a vector
   * stored at a heap address.
   */
  static const int NotFound;

 protected:
  /**
   * @see waitUntilAllBackgroundTasksHaveTerminated()
   */
  static int                                _NumberOfTriggeredTasks;

  /**
   * @see waitUntilAllBackgroundTasksHaveTerminated()
   */
  static tarch::multicore::BooleanSemaphore _heapSemaphore;

  /**
   * Each solver has an identifier/name. It is used for debug purposes only.
   */
  const std::string _identifier;

  const Type _type;

  /**
   * The number of state variables of the conservation or balance law.
   */
  const int _numberOfVariables;

  /**
   * The number of parameters, e.g, material parameters.
   */
  const int _numberOfParameters;

  /**
   * The number of nodal basis functions that are employed in each
   * coordinate direction.
   */
  const int _nodesPerCoordinateAxis;

  /**
   * The maximum extent a cell is allowed to have in each coordinate direction.
   *
   * \note This is an upper bound specified in the specification file.
   * This is not the actual maximum extent of a cell.
   */
  const double _maximumMeshSize;

  /**
   * The minimum extent in each coordinate direction at least one cell in the grid has.
   *
   * This value needs to be updated every time the grid has been changed.
   */
  double _minCellSize; // TODO(Dominic): Remove these cell size variables and put them into the solver subclasses.
  double _nextMinCellSize; // TODO(Dominic): Remove these cell size variables and put them into the solver subclasses.

  /**
   * The maximum extent in each coordinate direction at least one cell in the grid has.
   *
   * * This value needs to be updated every time the grid has been changed.
   */
  double _maxCellSize; // TODO(Dominic): Remove these cell size variables and put them into the solver subclasses.
  double _nextMaxCellSize; // TODO(Dominic): Remove these cell size variables and put them into the solver subclasses.

  /**
   * The time stepping mode of this solver.
   */
  const TimeStepping _timeStepping;

  /**
   * A profiler for this solver.
   */
  std::unique_ptr<profilers::Profiler> _profiler;

 public:
  Solver(const std::string& identifier, exahype::solvers::Solver::Type type,
         int numberOfVariables, int numberOfParameters,
         int nodesPerCoordinateAxis, double maximumMeshSize,
         exahype::solvers::Solver::TimeStepping timeStepping,
         std::unique_ptr<profilers::Profiler> profiler =
             std::unique_ptr<profilers::Profiler>(
                 new profilers::simple::NoOpProfiler("")));

  virtual ~Solver() { _profiler->writeToConfiguredOutput(); }

  // Disallow copy and assignment
  Solver(const Solver& other) = delete;
  Solver& operator=(const Solver& other) = delete;

  /**
   * Return a string representation for the type \p param.
   */
  static std::string toString(const exahype::solvers::Solver::Type& param);

  /**
   * Return a string representation for the time stepping mode \p param.
   */
  static std::string toString(const exahype::solvers::Solver::TimeStepping& param);

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

  virtual void updateNextMinCellSize(double minCellSize) {
    _nextMinCellSize = std::min( _nextMinCellSize, minCellSize );
  }

  virtual void updateNextMaxCellSize(double maxCellSize) {
    _nextMaxCellSize = std::max( _nextMaxCellSize, maxCellSize );
  }

  virtual double getNextMinCellSize() const {
    return _nextMinCellSize;
  }

  virtual double getNextMaxCellSize() const {
    return _nextMaxCellSize;
  }

  virtual double getMinCellSize() const {
    return _minCellSize;
  }

  virtual double getMaxCellSize() const {
    return _maxCellSize;
  }

  /**
   * Returns the identifier of this solver.
   */
  std::string getIdentifier() const;

  /**
   * Returns the type of this solver.
   */
  Type getType() const;

  /**
   * Returns the time stepping algorithm this solver is using.
   */
  TimeStepping getTimeStepping() const;

  /**
   * Returns the number of state variables.
   */
  int getNumberOfVariables() const;

  /**
   * Returns the number of parameters, e.g.,material constants etc.
   */
  int getNumberOfParameters() const;

  /**
   * If you use a higher order method, then this operation returns the
   * polynomial degree plus one. If you use a Finite Volume method, it
   * returns the number of cells within a patch per coordinate axis.
   */
  int getNodesPerCoordinateAxis() const;

  virtual std::string toString() const;

  virtual void toString(std::ostream& out) const;

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  virtual double getMinTimeStamp() const = 0;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  virtual double getMinTimeStepSize() const = 0;

  virtual void updateMinNextTimeStepSize(double value) = 0;

  virtual void initSolverTimeStepData(double value) = 0;

  /**
   * Copies the time stepping data from the global solver onto the patch's time
   * meta data.
   *
   * \param[in] element Index of the cell description in
   *                    the array at address \p cellDescriptionsIndex.
   */
  virtual void synchroniseTimeStepping(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * This routine is called if we perform
   * a time step update
   */
  virtual void startNewTimeStep() = 0;

  /**
   * This routine is called after the
   * mesh has been updated and we might have
   * to reinitialise some of the time
   * stamps or time step sizes before
   * we call startNewTimeStep().
   *
   * At the time of the call of this
   * function, a new next time step size
   * has already been computed.
   */
  virtual void reinitialiseTimeStepData() = 0;

  virtual double getMinNextTimeStepSize() const=0;

  /**
   * Run over all solvers and identify the minimal time stamp.
   */
  static double getMinSolverTimeStampOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal sum of minimal time stamp
   * plus the minimal time step size.
   *
   * The result is a lower bound of the minimum time stamp
   * that will be obtained in the following time step.
   */
  static double estimateMinNextSolverTimeStampOfAllSolvers();

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  static double getMinSolverTimeStepSizeOfAllSolvers();

  /**
   * Run over all solvers and identify the maximal time stamp.
   *
   * On the individual patches, we do use the min time stamp so
   * far, so the routine returns the maximum over all min solver
   * time stamps.
   */
  static double getMaxSolverTimeStampOfAllSolvers();

  static bool allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping scheme);

  static double getCoarsestMeshSizeOfAllSolvers();
  static double getFinestMaximumMeshSizeOfAllSolvers();

  /**
   * Run over all solvers and identify the maximum depth of adaptive
   * refinement employed.
   *
   * This number directly correlates with the number
   * of grid iterations to run for performing an erasing operation.
   */
  static int getMaxAdaptiveRefinementDepthOfAllSolvers();

  /**
   * Returns true if the index \p cellDescriptionsIndex
   * is a valid heap index.
   */
  virtual bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const = 0;

  /**
   * If an entry for this solver exists,
   * return the element index of the cell description
   * in the array at address \p cellDescriptionsIndex.
   * Otherwise and if \p cellDescriptionsIndex is an
   * invalid index, return Solver::NotFound.
   */
  virtual int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const = 0;

  /**
   * \see exahype::amr::computeSubcellPositionOfCellOrAncestor
   */
  virtual SubcellPosition computeSubcellPositionOfCellOrAncestor(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Modify a cell description in enter cell event.
   * This event should be used for single cell operations
   * like marking for refinement, erasing, augmenting,
   * or deaugmenting.
   *
   * Returns true if a performed action requires to
   * refine the mesh.
   *
   * \note We use this at the moment only
   * for refinement events. We can consider later
   * on to merge the time stepping functionality
   * (solution update, predictor comp.) into
   * this hook.
   */
  virtual bool enterCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /**
   * Refinement routine that should be used for
   * collective children-parent operations.
   *
   * Returns true if the grid cell associated with
   * the cell description can be erased due to
   * an erasing/deaugmenting process that was finished
   * for the cell description.
   * Returns false otherwise.
   *
   * \note We use this at the moment only
   * for refinement events. We can consider later
   * on to merge the time stepping functionality
   * (solution update, predictor comp.) into
   * this hook.
   */
  virtual bool leaveCell(
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
      exahype::Vertex* const coarseGridVertices,
      const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
      exahype::Cell& coarseGridCell,
      const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
      const int solverNumber) = 0;

  /////////////////////////////////////
  // CELL-LOCAL
  /////////////////////////////////////
  /**
   * Compute and return a new admissible time step
   * size for the cell description
   * \p element in the array at address
   * \p cellDescriptionsIndex.
   *
   * Then, update the time stamps and time step
   * sizes on the cell description accordingly.
   *
   * \return The new admissible time step size if the
   * cell description is of type Cell or
   * std::numeric_limits<double>::max().
   *
   * \note The update of the time stamps
   * and the time step sizes of the cell description
   * performed in this method can be revoked by the
   * solver's time stepping synchronisation mode.
   * If a time stepping mode like global or globalfixed
   * is chosen, then the changes made here are simply
   * overwritten. You might want to take a look
   * into method synchroniseTimeStepping().
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   * not copied once for each thread.
   *
   * \see startNewTimeStep(),
   *      synchroniseTimeStepping(int,int),
   *      synchroniseTimeStepping(CellDescription&)
   */
  virtual double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) = 0;

  /**
   * Impose initial conditions.
   *
   * \note Make sure to reset neighbour merge
   * helper variables in this method call.
   */
  virtual void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) = 0;

  /**
   * Update the solution of a cell description.
   *
   * \note Make sure to reset neighbour merge
   * helper variables in this method call.
   */
  virtual void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) = 0;

  /**
     * In this method, the solver can perform post-processing
     * operations, e.g., compression of the
     * degrees of freedoms.
     *
     * \param[in] element Index of the cell description
     *                    holding the data to send out in
     *                    the array at address \p cellDescriptionsIndex.
     *                    This is not the solver number.
     *
     * \see tryGetElement
     */
    virtual void preProcess(
        const int cellDescriptionsIndex,
        const int element) = 0;

  /**
   * In this method, the solver can perform post-processing
   * operations, e.g., compression of the
   * degrees of freedoms.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void postProcess(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
    * Prolongates face data from a parent ADERDGCellDescription to
    * \p cellDescription if the fine grid cell associated with
    * \p cellDescription is adjacent to a boundary of the
    * coarse grid cell associated with the parent cell description.
    *
    * Further zero out the face data of ancestors.
    *
    * \note This function assumes a top-down traversal of the grid and must thus
    * be called from the enterCell(...) or descend(...) mapping methods.
    *
    * \note Calling this method makes only sense if \p cellDescription is of type
    * Descendant or Ancestor.
    */
  virtual void prolongateDataAndPrepareDataRestriction(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   *
   *
   * \note This function assumes a bottom-up traversal of the grid and must thus
   * be called from the leaveCell(...) or ascend(...) mapping methods.
   */
  virtual void restrictData(
        const int cellDescriptionsIndex,
        const int element,
        const int parentCellDescriptionsIndex,
        const int parentElement,
        const tarch::la::Vector<DIMENSIONS,int>& subcellIndex) = 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   *
   * <h2>Temporary variables</h2>
   * See  exahype::mappings::Merging::prepareTemporaryVariables()
   * exahype::mappings::SolutionRecomputation::prepareTemporaryVariables()
   * for details on the size of the allocation of
   * the temporary variables.
   */
  virtual void mergeNeighbours(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2,
        double**                                  tempFaceUnknowns,
        double**                                  tempStateSizedVectors,
        double**                                  tempStateSizedSquareMatrices) = 0;

  /**
   * Take the cell descriptions \p element
   * from array at address \p cellDescriptionsIndex
   * and merge it with boundary data.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    at address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   *
   * <h2>Temporary variables</h2>
   * See  exahype::mappings::Merging::prepareTemporaryVariables()
   * exahype::mappings::SolutionRecomputation::prepareTemporaryVariables()
   * for details on the size of the allocation of
   * the temporary variables.
   */
  virtual void mergeWithBoundaryData(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
        double**                                  tempFaceUnknowns,
        double**                                  tempStateSizedVectors,
        double**                                  tempStateSizedSquareMatrices) =0;

  /**
   * After all merges with neighbour and boundary data have been performed for all the surrounding
   * faces of a cell (description), we prepare the neighbour merge helper variables
   * for the next neighbour merge in one of the following grid traversals.
   */
  virtual void prepareNextNeighbourMerging(
        const int cellDescriptionsIndex,
        const int element,
        exahype::Vertex* const fineGridVertices,
        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const = 0;

  #ifdef Parallel
  /**
   * Send solver data to neighbour rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * \see tryGetElement
   */
  virtual void sendDataToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to neighbour rank.
   */
  virtual void sendEmptyDataToNeighbour(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex with metadata.
   *
   * Currently, the neighbour metadata is only the neighbour
   * type as int \p neighbourTypeAsInt.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithNeighbourMetadata(
      const int neighbourTypeAsInt,
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Receive solver data from neighbour rank and write
   * it on the cell description \p element in
   * the cell descriptions array stored at \p
   * cellDescriptionsIndex.
   *
   * \note Peano only copies the calling mapping
   * for each thread. The solvers in the solver registry are
   *  not copied once for each thread.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   *
   * <h2>Temporary variables</h2>
   * See  exahype::mappings::Merging::prepareTemporaryVariables()
   * exahype::mappings::SolutionRecomputation::prepareTemporaryVariables()
   * for details on the size of the allocation of
   * the temporary variables.
   */
  virtual void mergeWithNeighbourData(
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
      const int                                    level) = 0;

  /**
   * Drop solver data from neighbour rank.
   */
  virtual void dropNeighbourData(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send solver data to master or worker rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to master or worker rank
   * due to fork or join.
   */
  virtual void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                    toRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from master or worker rank
   * that was sent out due to a fork or join. Wrote the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from master or worker rank
   * that was sent out due to a fork or join.
   */
  virtual void dropWorkerOrMasterDataDueToForkOrJoin(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  ///////////////////////////////////
  // WORKER->MASTER
  ///////////////////////////////////

  /**
   * Determine if the solver has to reduce data for the cell description
   * \p element at heap address \p cellDescriptionsIndex.
   * This is the case if we have adaptive mesh refinement
   * activated and a cell description located at the
   * master worker boundary is of type Ancestor (not EmptyAncestor).
   * In this case the worker will restrict face data up to the
   * master in this iteration.
   *
   * Note that this function must be called in
   * prepareSendToWorker(...) if you plug it into the Peano engine.
   *
   * Note that this function is not in control of determining when to reduce
   * time step data. This should be done outside of this function.
   */
  virtual bool hasToSendDataToMaster(
      const int cellDescriptionsIndex,
      const int element) = 0;

  /**
   * Send data to the master that is not
   * depending on a particular cell description.
   *
   * This operation might be used for the reduction of a global
   * minimum time step size over all MPI ranks.
   *
   * \note We always assume that
   * startNewTimeStep() has been already called on the
   * local solver instance. You thus
   * have to return the updated local time step size.
   *
   * \see startNewTimeStep(), mergeWithWorkerData()
   */
  virtual void sendDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from worker rank.
   */
  virtual void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send solver data to master rank. Read the data from
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToMaster(
      const int                                    masterRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to master rank.
   */
  virtual void sendEmptyDataToMaster(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from worker rank.
   * Write the data to the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithWorkerData(
      const int                                    workerRank,
      const int                                    workerTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from worker rank.
   */
  virtual void dropWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  ///////////////////////////////////
  // MASTER->WORKER
  ///////////////////////////////////

  /**
   * Send data to the worker that is not
   * depending on a particular cell description.
   *
   * This operation might be used for the synchronisation
   * of a global minimum time step size over all MPI ranks.
   *
   * \see startNewTimeStep(), mergeWithMasterData()
   */
  virtual void sendDataToWorker(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from master rank.
   */
  virtual void mergeWithMasterData(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send solver data to worker rank. Read the data from
   * the cell description \p element in the cell descriptions
   * vector stored at \p cellDescriptionsIndex.
   *
   * \param[in] element Index of the ADERDGCellDescription
   *                    holding the data to send out in
   *                    the heap vector at \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void sendDataToWorker(
      const int                                    workerRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Send empty solver data to worker rank.
   */
  virtual void sendEmptyDataToWorker(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Merge with solver data from master rank
   * that was sent out due to a fork or join. Write the data to
   * the cell description \p element in
   * the cell descriptions vector stored at \p
   * cellDescriptionsIndex.
   *
   * \param[in] element Index of the cell description
   *                    holding the data to send out in
   *                    the array with address \p cellDescriptionsIndex.
   *                    This is not the solver number.
   */
  virtual void mergeWithMasterData(
      const int                                    masterRank,
      const int                                    masterTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;

  /**
   * Drop solver data from master rank.
   */
  virtual void dropMasterData(
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) = 0;
  #endif
};

#endif
