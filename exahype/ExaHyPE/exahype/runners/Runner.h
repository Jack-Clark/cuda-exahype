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
 
#ifndef _EXAHYPE_RUNNERS_RUNNER_H_
#define _EXAHYPE_RUNNERS_RUNNER_H_

#include "exahype/Parser.h"
#include "tarch/logging/Log.h"

#include "exahype/State.h"

namespace exahype {
namespace runners {
class Runner;
}
namespace repositories {
class Repository;
}
}

/**
 * Runner
 *
 */
class exahype::runners::Runner {
 private:
  static tarch::logging::Log _log;

  const exahype::Parser& _parser;

  /**
   * Setup the oracles for the shared memory parallelisation. Different
   * oracles can be employed:
   *
   * - If no autotuning is used and no valid properties file is provided and
   *   the code is compiled with -DPerformanceAnalysis, we use the grain size
   *   sampling
   * - If no autotuning is used and no valid properties file is provided, we
   *   use the default oracle coming along with the Peano kernel
   * - If autotuning is enabled and no valid properties file is provided, we
   *
   *
   * <h2>Invocation sequence</h2>
   *
   * It is important that we init the shared memory environment after we have
   * created the repository. See Orace::loadStatistics().
   */
  void initSharedMemoryConfiguration();

  /**
   * The shared memory environment has to be set up before we create the
   * repository.
   */
  void initDistributedMemoryConfiguration();
  void shutdownSharedMemoryConfiguration();
  void shutdownDistributedMemoryConfiguration();

  int runAsMaster(exahype::repositories::Repository& repository);

#ifdef Parallel
  int runAsWorker(exahype::repositories::Repository& repository);

  /**
   * If the master node calls runGlobalStep() on the repository, all MPI
   * ranks besides rank 0 invoke this operation no matter whether they are
   * idle or not. Hence, please call this operation manually within
   * runAsMaster() if you require the master to execute the same function
   * as well.
   */
  void runGlobalStep();
#endif

  /**
   * Reset all time stamps to zero. Runs through the solver registry only,
   * i.e. no grid traversal is required.
   */
  void initSolverTimeStepData();

  void validateInitialSolverTimeStepData(const bool fuseADERDGPhases) const;

  /**
   * Initialise the data compression (or switch it off if we don't need it).
   * The routine is called for each and every rank.
   */
  void initDataCompression();


  /**
   * Print minimum of current solver time stamps and time step sizes.
   *
   * The solver time stamp and step sizes are computed as
   * minimum over the patch quantities.
   *
   * \param numberOfStepsRanSinceLastCall The number of time step done since
   *          the last call. If you pass -1, then you've done only preparatory
   *          time steps so far.
   *
   * \note A value \p numberOfStepsRanSinceLastCall greater than 1 is only interesting
   *       for fixed time stepping runs.
   */
  void printTimeStepInfo(int numberOfStepsRanSinceLastCall, const exahype::repositories::Repository& repository);


  /**
   * Do one time step where all phases are actually fused into one traversal
   *
   * @param numberOfStepsToRun Number of steps to run. If you hand in 0, then
   *           it runs one time step plus does a plot.
   */
  void runOneTimeStampWithFusedAlgorithmicSteps(exahype::repositories::Repository& repository, int numberOfStepsToRun, bool exchangeBoundaryData);

  /**
   * If the repository state indicates that a time update did
   * harm an ADER-DG solver's stability condition, we rerun the
   * prediction phase (for all ADER-DG solvers for now).
   *
   * \note Currently this function only makes sense for single
   * ADER-DG solvers and global time stepping.
   */
  void recomputePredictorIfNecessary(
      exahype::repositories::Repository& repository);

  /**
   * Run the three adapters necessary for updating the
   * limiter domain.
   *
   * \param[in] reinitialiseTimeStepData This flag must be set when we initialise the limiter domain.
   *
   */
  void updateLimiterDomain(exahype::repositories::Repository& repository);

  /**
   * Do one time step but actually use a couple of iterations to do so.
   *
   *
   * @param plot      Do plot in the after the corrector has been applied
   */
  void runOneTimeStampWithThreeSeparateAlgorithmicSteps(
      exahype::repositories::Repository& repository, bool plot);

  void validateSolverTimeStepDataForThreeAlgorithmicPhases(const bool fuseADERDGPhases) const;

  /**
   * Sets up the geometry, hands it over to a new instance of the repository
   * and returns the repository.
   */
  exahype::repositories::Repository* createRepository() const;

  /**
   * Constructs the initial computational grid
   *
   * The grid generation is an iterative process as it may happen that a grid
   * is refined by a cell initialisation. Cell-based refinements however never
   * are realised in the exactly same iteration, so we always have to wait yet
   * another iteration.
   *
   * This solves immediately another problem: ExaHyPE's adjacency management
   * requires us to run once more over the grid at least once its completely
   * built up to get all the adjacency information correct. If we rely on the
   * grid to become stationary, this is always the case - as long as additional
   * vertices are added, the grid remains instationary. When we've added the
   * last grid entities and run the adapter once again, then it becomes
   * stationary.
   *
   * For the parallel case, I've changed from stationary into balanced which is
   * a slight generalisation. See Peano guidebook.
   */
  void createGrid(exahype::repositories::Repository& repository);

  /**
   * Run through all the solvers and identify the coarsest grid level in the tree
   * that will be populated by a solver.
   */
  int getCoarsestGridLevelOfAllSolvers() const;
 public:
  explicit Runner(const Parser& parser);
  virtual ~Runner();

  // Disallow copy and assignment
  Runner(const Runner& other) = delete;
  Runner& operator=(const Runner& other) = delete;

  /**
   * Run
   */
  int run();
};

#endif
