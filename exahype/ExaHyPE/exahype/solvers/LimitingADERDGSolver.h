/*
 * LimitingADERDGSolver.h
 *
 *  Created on: 26 Oct 2016
 *      Author: dominic
 */

#ifndef LIMITEDADERDGSOLVER_H_
#define LIMITEDADERDGSOLVER_H_

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/FiniteVolumesSolver.h"

#include "tarch/multicore/BooleanSemaphore.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/mappings/Merging.h"
#include "exahype/mappings/Prediction.h"
#include "exahype/mappings/LimiterStatusSpreading.h"
#include "exahype/mappings/Reinitialisation.h"
#include "exahype/mappings/SolutionRecomputation.h"

#include "exahype/plotters/Plotter.h"

#include "exahype/profilers/simple/NoOpProfiler.h"

namespace exahype {
namespace solvers {

class LimitingADERDGSolver;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::LimitingADERDGSolver : public exahype::solvers::Solver {

private:
  typedef exahype::records::ADERDGCellDescription SolverPatch;
  typedef peano::heap::PlainHeap<SolverPatch> SolverHeap;

  typedef exahype::records::FiniteVolumesCellDescription LimiterPatch;
  typedef peano::heap::PlainHeap<LimiterPatch> LimiterHeap;


  /**
   * Log device.
   */
  static tarch::logging::Log _log;

  /**
   * The ADERDG solver.
   */
  std::unique_ptr<exahype::solvers::ADERDGSolver> _solver;

  /**
   * The finite volumes solver used for the a posteriori subcell limiting.
   */
  std::unique_ptr<exahype::solvers::FiniteVolumesSolver> _limiter;

  /**
   * A flag indicating that the limiter domain has changed.
   * This might be the case if either a cell has been
   * newly marked as troubled or healed.
   */
  bool _limiterDomainHasChanged;

  /**
   * The limiterDomainHasChanged for the next
   * iteration.
   */
  bool _nextLimiterDomainHasChanged;

  /**
   * The maximum relaxation parameter
   * used for the discrete maximum principle.
   */
  double _DMPMaximumRelaxationParameter;

  /**
   * The difference scaling
   * used for the discrete maximum principle.
   */
  double _DMPDifferenceScaling;

  /**
   * TODO(Dominc): Remove after docu is recycled.
   *
   * This operation sets the solutions' minimum and maximum value on a cell.
   * The routine is to be invoked after the code has determined the new minimum
   * and maximum value within a cell. In turn, it evaluates whether the new
   * minimum and maximum value have decreased or grown, respectively.
   *
   * If the new min/max values indicate that the new solution comprises
   * oscillations, the routine returns false. This is an indicator that the
   * solution should be limited.
   *
   * If the new min/max values fit, the routine returns true.
   *
   * <h2>Implementation</h2>
   * We hold the min/max information exclusively on the faces. The first thing
   * the routine does is to project the min/max values into the cell. For this
   * it evaluates the 2d faces. The projected value then is compared to the
   * arguments. Once the results of the operation is determined, the routine
   * writes the new arguments onto the 2d face entries. This, on the one hand,
   * stores the data for the subsequent time step, but it also propagates the
   * min/max information into the face-connected neighbours.
   *
   * @param  min          New minimum values within the cell. Array of length
   *                      _numberOfUnknowns.
   * @param  max          New maximum values within the cell
   * @param  solverIndex  Number of the solver within the cell. Please ensure
   *                      that solverIndex refers to an ADER-DG solver.
   * @return True if the new min and max values fit into the restricted min
   *   max solutions. Return false if we seem to run into oscillations.
   */
  //  void setSolutionMinMax(double* min, double* max) const;

  /**
   * Merge the solution min and max values on a face between two cell
   * descriptions. Signature is similar to that of the solver of a Riemann problem.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch& pLeft,
      SolverPatch& pRight,
      const int faceIndexLeft,
      const int faceIndexRight) const;

  /**
   * Determine a new merged limiter status based on the neighbours merged
   * limiter status.
   *
   * This method is used during the limiter status spreading.
   * It assumes that the merged limiter statuses are initialised with
   * Ok or Troubled before the spreading.
   * It further assumes that a unique value is written again to
   * the merged limiter status fields after the limiter status of all neighbours
   * of a solver patch have been merged with the limiter status of the patch.
   * (see exahype::solvers::LimitingADERDG::updateLimiterStatus).
   *
   * Determining the merged limiter status:
   * | Status     | Neighbour's Status | Merged Status |
   * --------------------------------------------------|
   * | O          | O/NNT              | O
   * | ~          | T                  | NT
   * | ~          | NT                 | NNT
   * | NNT        | T                  | NT
   */
  void mergeWithNeighbourLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus) const;

  void mergeWithNeighbourLimiterStatus(
      SolverPatch& solverPatch,
      const int faceIndex,
      const SolverPatch::LimiterStatus& neighbourLimiterStatus,
      const SolverPatch::LimiterStatus& neighbourOfNeighbourLimiterStatus) const;

  /**
   * Determine a unified limiter status of a solver patch after a limiter status spreading
   * iteration.
   *
   * <h2>Determining the unified value</h2>
   * If all of the merged limiter status fields
   * are set to Troubled, the limiter status is changed to Troubled.
   * (There is either all or none of the statuses set to Troubled.)
   *
   * Otherwise and if at least one of the merged statuses is set to NeighbourOfTroubledCell,
   * the status is set to NeighbourOfTroubledCell.
   *
   * Otherwise and if at least one of the merged statuses is set to NeighbourIsNeighbourOfTroubledCell,
   * the status is set to NeighbourIsNeighbourOfTroubledCell.
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   */
  exahype::solvers::LimitingADERDGSolver::SolverPatch::LimiterStatus determineLimiterStatus(SolverPatch& solverPatch) const;

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   *
   * Compute the new min and max at the same time.
   */
  bool evaluateDiscreteMaximumPrincipleAndDetermineMinAndMax(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver contains unphysical oscillations (false)
   * or not (true).
   */
  bool evaluateDiscreteMaximumPrinciple(SolverPatch& solverPatch);

  /**
   * Checks if the updated solution
   * of the ADER-DG solver is
   * a physically admissible one (true).
   */
  bool evaluatePhysicalAdmissibilityCriterion(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the solution (per variable)
   * and makes them accessible per face.
   */
  void determineSolverMinAndMax(SolverPatch& solverPatch);

  /**
   * Computes the cell-local minimum and maximum
   * of the limiter's solution (per variable)
   * and makes them accessible per face.
   *
   * This method is used for troubled cells that
   * do not hold a valid ADER-DG solution,
   * as well as their neighbours.
   */
  void determineLimiterMinAndMax(SolverPatch& solverPatch,LimiterPatch& limiterPatch);

  /**
   * Updates the merged limiter status based on the cell-local ADER-DG solution
   * values,
   */
  void updateMergedLimiterStatusAfterSolutionUpdate(SolverPatch& solverPatch,const bool isTroubled);

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
   * Send the solution minimum and maximum values per variable
   * and further the merged limiter status of the solver patch
   * \p element in heap array \p cellDescriptionsIndex to the
   * neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void sendMinAndMaxToNeighbour(
      const int                                    toRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Receive the solution minimum and maximum values per variable
   * and further the merged limiter status for the solver patch
   * \p element in heap array \p cellDescriptionsIndex from the
   * neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for a description of the parameters.
   */
  void mergeWithNeighbourMinAndMax(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Single-sided variant of mergeSolutionMinMaxOnFace() that is required
   * for MPI where min and max value are explicitly exchanged through messages.
   */
  void mergeSolutionMinMaxOnFace(
      SolverPatch&  cellDescription,
      const int     faceIndex,
      const double* const min, const double* const  max) const;

//  /**
//   * Send the limiter status
//   * of the solver patch \p element in heap array
//   * \p cellDescriptionsIndex to the respective neighbour.
//   *
//   * \see exahype::solvers::Solver::sendDataToNeighbour
//   * for a description of the parameters.
//   *
//   * <h2>Heap</h2>
//   * We currently use the DoubleHeap to
//   * communicate the enum.
//   * While it would make sense to use the
//   * MetadataHeap for this purpose,
//   * there is a possibility of intermixing
//   * the communication of the limiter status
//   * with the communication of the solver
//   * metadata.
//   */
//  void sendMergedLimiterStatusToNeighbour(
//      const int                                     toRank,
//      const int                                     cellDescriptionsIndex,
//      const int                                     element,
//      const tarch::la::Vector<DIMENSIONS, int>&     src,
//      const tarch::la::Vector<DIMENSIONS, int>&     dest,
//      const tarch::la::Vector<DIMENSIONS, double>&  x,
//      const int                                     level) const;
#endif

public:
  static bool limiterDomainOfOneSolverHasChanged();

  /**
   * TODO(Dominic): Docu
   *
   *
   * <h2>Discrete maximum principle</h2>
   * By default this constructor initialises the maximum relaxation
   * parameter to the value to \f$ \delta_0 = 1\cdot 10^{-4} \f$
   * and the difference scaling parameter to \f$ \epsilon = 1\cdot 10^{-3} \f$.
   * See Dumbser et al., 2014. doi:10.1016/j.jcp.2014.08.009 for more details on
   * the notation.
   */
  LimitingADERDGSolver(
      const std::string& identifier,
      std::unique_ptr<exahype::solvers::ADERDGSolver> solver,
      std::unique_ptr<exahype::solvers::FiniteVolumesSolver> limiter,
      const double DMPRelaxationParameter=1e-4,
      const double DMPDifferenceScaling=1e-3
      );

  virtual ~LimitingADERDGSolver() {
    _solver.reset();
    _limiter.reset();
  }

  // Disallow copy and assignment
  LimitingADERDGSolver(const ADERDGSolver& other) = delete;
  LimitingADERDGSolver& operator=(const ADERDGSolver& other) = delete;

  /*
   * A time stamp minimised over all the ADERDG and FV solver
   * patches.
   */
  double getMinTimeStamp() const override;

  /**
   * Run over all solvers and identify the minimal time step size.
   */
  double getMinTimeStepSize() const override;

  double getMinNextTimeStepSize() const override;

  void updateMinNextTimeStepSize(double value) override;

  void initSolverTimeStepData(double value) override;

  void synchroniseTimeStepping(
          const int cellDescriptionsIndex,
          const int element) override;

  void reconstructStandardTimeSteppingData() {
    _solver->reconstructStandardTimeSteppingData();
    _limiter->setMinTimeStamp(_solver->getMinCorrectorTimeStamp());
    _limiter->setMinTimeStepSize(_solver->getMinCorrectorTimeStepSize());
  }

  void startNewTimeStep() override;

  bool getLimiterDomainHasChanged() {
    return _limiterDomainHasChanged;
  }

  bool getNextLimiterDomainHasChanged() {
    return _nextLimiterDomainHasChanged;
  }

  void updateNextLimiterDomainHasChanged(bool state) {
    _nextLimiterDomainHasChanged |= state;
  }

  /**
   * Roll back the time step data to the
   * ones of the previous time step.
   */
  void rollbackToPreviousTimeStep();

  void reconstructStandardTimeSteppingDataAfterRollback();

  void reinitialiseTimeStepData() override;

  void updateNextMinCellSize(double minCellSize) override;

  void updateNextMaxCellSize(double maxCellSize) override;

  double getNextMinCellSize() const override;

  double getNextMaxCellSize() const override;

  double getMinCellSize() const override;

  double getMaxCellSize() const override;

  bool isValidCellDescriptionIndex(
      const int cellDescriptionsIndex) const override ;

  /**
   * Returns the index of the solver patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const override {
    return _solver->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Returns the index of the limiter patch registered for the solver with
   * index \p solverNumber in exahype::solvers::RegisteredSolvers.
   * If no patch is registered for the solver, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElement(
      const int cellDescriptionsIndex,
      const int solverNumber) const {
    return _limiter->tryGetElement(cellDescriptionsIndex,solverNumber);
  }

  /**
   * Returns the index of the limiter patch corresponding to
   * the solver patch with index \p solverElement.
   * Both patches link to the same ::LimitingADERDGSolver.
   * If no limiter patch is found, this function
   * returns exahype::solvers::Solver::NotFound.
   */
  int tryGetLimiterElementFromSolverElement(
      const int cellDescriptionsIndex,
      const int solverElement) const {
    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,solverElement);
    return _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
  }

  /**
    * \see exahype::amr::computeSubcellPositionOfCellOrAncestor
    */
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
  //////////////////////////////////
  double startNewTimeStep(
      const int cellDescriptionsIndex,
      const int element,
      double*   tempEigenvalues) override;


  void reconstructStandardTimeSteppingData(const int cellDescriptionsIndex,int element) const {
    _solver->reconstructStandardTimeSteppingData(cellDescriptionsIndex,element);

    SolverPatch& solverPatch = _solver->getCellDescription(cellDescriptionsIndex,element);
    const int limiterElement = _limiter->tryGetElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
    if (limiterElement!=exahype::solvers::Solver::NotFound) {
      LimiterPatch& limiterPatch = _limiter->getCellDescription(cellDescriptionsIndex,limiterElement);
      limiterPatch.setTimeStamp(solverPatch.getCorrectorTimeStamp());
      limiterPatch.setTimeStepSize(solverPatch.getCorrectorTimeStepSize());
    }
  }

 /**
   * Rollback to the previous time step, i.e,
   * overwrite the time step size and time stamp
   * fields of the solver and limiter patches
   * by the values used in the previous iteration.
   */
  void rollbackToPreviousTimeStep(
      const int cellDescriptionsIndex,
      const int solverElement);

  /**
   * Similar to reconstructStandardTimeSteppingData for roll backs
   */
  void reconstructStandardTimeSteppingDataAfterRollback(
      const int cellDescriptionsIndex,
      const int solverElement) const;

  /**
   * TODO(Dominic): I need the whole limiter recomputation
   * procedure also for the initial conditions.
   *
   * This includes computing, sending, and merging
   * of the min/max values.
   */
  void setInitialConditions(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * Adds a limiter patch to the initially troubled cells, their
   * neighbours, and their neighbour's neighbours.
   * Further imposes finite volumes boundary conditions onto the limiter
   * solution in cells with limiter status Troubled and NeighbourIsTroubled and
   * then projects this solution onto the DG solver space for those cells.
   * Finally, projects the ADER-DG initial conditions onto
   * the FV limiter space for cells with limiter status
   * NeighbourIsNeighbourOfTroubledCell.
   */
  void initialiseLimiter(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const;

  /**
   * This method assumes the ADERDG solver's cell-local limiter status has
   * already been determined.
   *
   * \see determineLimiterStatusAfterLimiterStatusSpreading(...)
   */
  void updateSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * Determine the new cell-local min max values.
   *
   * Must be invoked after ::determineLimiterStatusAfterSolutionUpdate.
   */
  void determineMinAndMax(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Evaluates a discrete maximum principle (DMP) and
   * the physical admissibility detection (PAD) criterion for
   * the solution values stored for a solver patch.
   *
   * This method then invokes
   * ::determinMergedLimiterStatusAfterSolutionUpdate(SolverPatch&,const bool)
   * with the result of these checks.
   */
  bool updateMergedLimiterStatusAndMinAndMaxAfterSolutionUpdate(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Similar to ::determineMergedLimiterStatusAfterSolutionUpdate(const int,const int)
   * but does not evaluate the discrete maximum principle.
   */
  bool updateMergedLimiterStatusAndMinAndMaxAfterSetInitialConditions(
      const int cellDescriptionsIndex,
      const int element);

  /**
   * Update the merged limiter status based on
   * the cell-local solution values.
   *
   * \param[in] isTroubled A bool indicating if the patch's solution is (still) troubled
   *
   * \return True if the limiter domain changes in the cell, i.e.,
   * if a patch with status Ok,NeighbourIsTroubled, or NeighbourOfNeighbourIsTroubled
   * changes its status to Troubled. Or, if a patch with status Troubled changes its
   * status to Ok.
   */
  bool determineMergedLimiterStatusAfterSolutionUpdate(
      SolverPatch& solverPatch,const bool isTroubled) const;

  /**
   * Update the merged limiter status with a unified value based on
   * the merged limiter statuses per face of the cell.
   *
   * \see ::determineLimiterStatus
   */
  void updateMergedLimiterStatus(
      const int cellDescriptionsIndex,const int element) const;

  /**
   * Update the limiter status with a unified value based on
   * the merged limiter statuses per face of the cell.
   *
   * \see ::determineLimiterStatus
   */
  void updateLimiterStatus(
      const int cellDescriptionsIndex,const int element) const;

  /**
   * Reinitialises cells that have been subject to a limiter status change.
   * This method is invoked (during and??) after the limiter status spreading.
   *
   * The method has to take into account which solution, the solver's
   * or the limiter's, was populated with valid solution values
   * in the last iteration. The action of this method is
   * thus based on the new and old limiter status.
   *
   * We perform the following actions based on the
   * old and new limiter status:
   *
   * | New Status | Old Status | Action                                                                                                                                        |
   * ------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
   * | O          | O          | Do nothing.                                                                                                                                   |
   * | ~          | T/NT/NNT   | Remove the limiter patch.                                                                                                                     |
   * | T/NT/NNT   | T/NT       | Roll back the limiter solution.                                                                                                               |
   * | ~          | O/NNT      | Roll back the solver solution. Initialise limiter patch if necessary. Project (old, valid) solver solution onto the limiter's solution space. |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * TODO(Dominic)
   * Adapters:
   * LimitingADERDGSolver LimiterStatusSpreading
   * LimitingADERDGSolver Reinitialisation
   * LimitingADERDGSolver SolutionRecomputation
   */
  void reinitialiseSolvers(
      const int cellDescriptionsIndex,
      const int element,
      exahype::Cell& fineGridCell,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const;

  /**
   * Recompute the solution in cells that have been subject to a limiter status change
   * This method is invoked after the solver reinitialisation
   * (see exahype::solvers::LimitingADERDGSolver::reinitialiseSolvers).
   *
   * It evolves the solution of the solver and limiter in the reinitialised cells to the
   * correct time stamp.
   *
   * We perform the following actions based on the
   * new limiter status:
   *
   * |New Status | Action                                                                                                                                      |
   * ----------------------------------------------------------------------------------------------------------------------------------------------------------|
   * |O          | Do nothing. Solver solution has been evolved correctly before.                                                                              |
   * |T/NT       | Evolve FV solver project result onto the ADER-DG space.                                                                                     |
   * |NNT        | Evolve solver and project its solution onto the limiter solution space. (We had to do a rollback beforehand in the reinitialisation phase.) |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * We do not overwrite the old limiter status set in this method.
   * We compute the new limiter status based on the merged limiter statuses associated
   * with the faces.
   *
   * <h2>Overlapping status spreading and reinitialisation with solution reconputation</h2>
   * We can recompute the new solution in cells with status Troubled after one iteration
   * since old solution values from direct neighbours are available then.
   *
   * We can recompute the
   *
   * TODO(Dominic)
   * Adapters:
   * LimitingADERDGSolver LimiterStatusSpreading
   * LimitingADERDGSolver Reinitialisation
   * LimitingADERDGSolver SolutionRecomputation
   */
  void recomputeSolution(
      const int cellDescriptionsIndex,
      const int element,
      double** tempStateSizedArrays,
      double** tempUnknowns,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator);

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
  // NEIGHBOUR
  ///////////////////////////////////
  /**
   *
   */
  void mergeLimiterStatusOfNeighbours(
        const int                                 cellDescriptionsIndex1,
        const int                                 element1,
        const int                                 cellDescriptionsIndex2,
        const int                                 element2,
        const tarch::la::Vector<DIMENSIONS, int>& pos1,
        const tarch::la::Vector<DIMENSIONS, int>& pos2);

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

  /**
   * Merge solver boundary data (and other values) of two adjacent
   * cells based on their limiter status.
   *
   * The solver involved in the neighbour merge
   * is selected according to the following scheme:
   *
   * | Status 1 | Status 2 | Solver to Merge
   * ---------------------------------------
   * | O        | O        | ADER-DG       |
   * | O        | NNT      | ADER-DG       |// O|NNT x O|NNT
   * | NNT      | O        | ADER-DG       |
   * | NNT      | NNT      | ADER-DG       |
   *
   * | NNT      | NT       | FV            |
   * | NT       | NNT      | FV            | // NT&NNT | N&NNT
   *
   * | NT       | NT       | FV            |
   * | NT       | T        | FV            |
   * | T        | NT       | FV            | // T|NT x T|NT
   * | T        | T        | FV            |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with neighbour data
   * in the recomputation phase.
   *
   * TODO(Dominic): Remove limiterstatus1 and limiterStatus2 argument.
   * They depend on the isRecomputation value
   */
  void mergeNeighboursBasedOnLimiterStatus(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      const bool                                isRecomputation,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices);

  /**
   * Merges only the min max of two neighbours.
   *
   * This method is used to detect cells that are
   * troubled after the imposition of initial conditions.
   */
  void mergeSolutionMinMax(
      const int                                 cellDescriptionsIndex1,
      const int                                 element1,
      const int                                 cellDescriptionsIndex2,
      const int                                 element2,
      const tarch::la::Vector<DIMENSIONS, int>& pos1,
      const tarch::la::Vector<DIMENSIONS, int>& pos2,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices);

  /**
   * Depending on the limiter status, we impose boundary conditions
   * onto the solution of the solver or of the limiter.
   */
  void mergeWithBoundaryData(
      const int                                 cellDescriptionsIndex,
      const int                                 element,
      const tarch::la::Vector<DIMENSIONS, int>& posCell,
      const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
      double**                                  tempFaceUnknowns,
      double**                                  tempStateSizedVectors,
      double**                                  tempStateSizedSquareMatrices) override;

  /**
   * Merge solver boundary data (and other values) of a
   * cell with the boundary conditions based on the cell's
   * limiter status.
   *
   * The solver involved in the merge
   * is selected according to the following scheme:
   *
   * | Status   | Solver to Merge |
   * ------------------------------
   * | O        | ADER-DG         |
   * | NNT      | ADER-DG         |
   *
   * | NT       | FV              |
   * | T        | FV              |
   *
   * Legend: O: Ok, T: Troubled, NT: NeighbourIsTroubledCell, NNT: NeighbourIsNeighbourOfTroubledCell
   *
   * <h2>Solution Recomputation</h2>
   * If we perform a solution recomputation, we do not need to perform
   * a solution update in the non-troubled solver patches and the patches
   * with status neighbourIsNeighbourOfTroubledCell. Instead we
   * can simply reuse the already computed update to advance in
   * time to the desired time stamp.
   *
   * We thus do not need to merge these patches with boundary data
   * in the recomputation phase.
   *
   * \param[in] isRecomputation Flag indicating if this merge is part of a solution recomputation phase.
   */
  void mergeWithBoundaryDataBasedOnLimiterStatus(
        const int                                 cellDescriptionsIndex,
        const int                                 element,
        const SolverPatch::LimiterStatus&         limiterStatus,
        const tarch::la::Vector<DIMENSIONS, int>& posCell,
        const tarch::la::Vector<DIMENSIONS, int>& posBoundary,
        const bool                                isRecomputation,
        double**                                  tempFaceUnknowns,
        double**                                  tempStateSizedVectors,
        double**                                  tempStateSizedSquareMatrices);

  void prepareNextNeighbourMerging(
      const int cellDescriptionsIndex,const int element,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) const override;

#ifdef Parallel
  ///////////////////////////////////
  // NEIGHBOUR
  ///////////////////////////////////
  void mergeWithNeighbourMetadata(
        const int neighbourTypeAsInt,
        const int cellDescriptionsIndex,
        const int element) override;

  void sendDataToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  /**
   * Send data or empty data to the neighbour data based
   * on the limiter status.
   *
   * \param[in] isRecomputation Indicates if this called within a solution recomputation
   *                            process.
   * \param[in] limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                            which either make use of the unified face-wise limiter status (isRecomputation)
   *                            or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>SolutionRecomputation</h2>
   * In case this method is called within a solution recomputation,
   * the ADER-DG solver does only send empty messages to the neighbour.
   * Otherwise it merges the received data and adds it to the update.
   *
   * \note This method assumes that there has been a unified face-wise limiter status value
   * determined and written back to the faces a-priori.
   *
   * <h2>Possible optimisations</h2>
   * Depending on isRecomputation we do not need to send both, solver and limiter,
   * data for patches with status NeighbourIsNeighbourOfTroubledCell and NeighbourOfTroubledCell.
   */
  void sendDataToNeighbourBasedOnLimiterStatus(
        const int                                    toRank,
        const int                                    cellDescriptionsIndex,
        const int                                    element,
        const tarch::la::Vector<DIMENSIONS, int>&    src,
        const tarch::la::Vector<DIMENSIONS, int>&    dest,
        const bool                                   isRecomputation,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const int                                    level) const;

  void sendEmptyDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

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
   * Merge or drop received neighbour data based
   * on the limiter status.
   *
   * \param isRecomputation Indicates if this called within a solution recomputation
   *                        process.
   * \param limiterStatus   This method is used in two different contexts (see \p isRecomputation)
   *                        which either make use of the unified face-wise limiter status (isRecomputation)
   *                        or make use of the cell-wise limiter status (!isRecomputation).
   *
   * <h2>SolutionRecomputation</h2>
   * In case this method is called within a solution recomputation,
   * the ADER-DG solver drops the received boundary data.
   * Otherwise it merges the received data and adds it to the update.
   *
   *  \note This method assumes that there has been a unified face-wise limiter status value
   *  determined and written back to the faces.
   */
  void mergeWithNeighbourDataBasedOnLimiterStatus(
      const int                                    fromRank,
      const int                                    neighbourTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const bool                                   isRecomputation,
      double**                                     tempFaceUnknowns,
      double**                                     tempStateSizedVectors,
      double**                                     tempStateSizedSquareMatrices,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  void dropNeighbourData(
      const int                                     fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  ///////////////////////////////////////
  // NEIGHBOUR - Limiter status spreading
  ///////////////////////////////////////

  /**
   * Receive and merge the merged limiter status sent
   * by the neighbour at position \p src - \p dest.
   *
   * see exahype::solvers::Solver::mergeWithNeighbourData
   * for details on the parameters.
   */
  void mergeWithNeighbourMergedLimiterStatus(
      const int                                    fromRank,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   * Drop the merged limiter status sent
   * by the neighbour at position \p src - \p dest.
   *
   * \see exahype::solvers::Solver::dropNeighbourData
   * for details on the parameters.
   */
  void dropNeighbourMergedLimiterStatus(
      const int                                    fromRank,
      const tarch::la::Vector<DIMENSIONS, int>&    src,
      const tarch::la::Vector<DIMENSIONS, int>&    dest,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) const;

  /**
   *  Send the merged limiter status at position \p dest- \p src to
   *  the corresponding neighbour.
   *
   * \see exahype::solvers::Solver::sendDataToNeighbour
   * for details on the parameters.
   */
  void sendMergedLimiterStatusToNeighbour(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  /**
   *  Send an empty array instead of
   *  the merged limiter status at position \p dest- \p src to
   *  the corresponding neighbour.
   *
   *  \see exahype::solvers::Solver::sendEmptyDataToNeighbour
   * for details on the parameters.
   */
  void sendEmptyDataInsteadOfMergedLimiterStatusToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  ///////////////////////////////////////
  // NEIGHBOUR - Solution Recomputation
  ///////////////////////////////////////
  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void sendEmptySolverAndLimiterDataToNeighbour(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, int>&     src,
      const tarch::la::Vector<DIMENSIONS, int>&     dest,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) const;

  /**
   * TODO(Dominic):
   * Add more docu.
   *
   * We do not send the minimum and maximum solution values
   * during the solution recomputation process.
   */
  void dropNeighbourSolverAndLimiterData(
        const int                                     fromRank,
        const tarch::la::Vector<DIMENSIONS, int>&     src,
        const tarch::la::Vector<DIMENSIONS, int>&     dest,
        const tarch::la::Vector<DIMENSIONS, double>&  x,
        const int                                     level) const;

  /////////////////////////////////////
  // FORK OR JOIN
  /////////////////////////////////////
  void sendDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorkerOrMasterDueToForkOrJoin(
      const int                                     toRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerOrMasterDataDueToForkOrJoin(
      const int                                     fromRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

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
      const int                                    masterRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void sendDataToMaster(
      const int                                     masterRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToMaster(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithWorkerData(
      const int                                    workerRank,
      const int                                    workerTypeAsInt,
      const int                                    cellDescriptionsIndex,
      const int                                    element,
      const tarch::la::Vector<DIMENSIONS, double>& x,
      const int                                    level) override;

  void dropWorkerData(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

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
      const int                                     workerRank,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void sendEmptyDataToWorker(
      const int                                     workerRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void mergeWithMasterData(
      const int                                     masterRank,
      const int                                     masterTypeAsInt,
      const int                                     cellDescriptionsIndex,
      const int                                     element,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;

  void dropMasterData(
      const int                                     masterRank,
      const tarch::la::Vector<DIMENSIONS, double>&  x,
      const int                                     level) override;
#endif

  std::string toString() const override;

  void toString (std::ostream& out) const override;

  const std::unique_ptr<exahype::solvers::FiniteVolumesSolver>&
  getLimiter () const {
    return _limiter;
  }

  const std::unique_ptr<exahype::solvers::ADERDGSolver>&
  getSolver () const {
    return _solver;
  }
};


#endif /* LIMITEDADERDGSOLVER_H_ */
