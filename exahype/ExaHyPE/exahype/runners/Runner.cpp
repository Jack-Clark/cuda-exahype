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

#include "exahype/runners/Runner.h"

#include "../../../Peano/mpibalancing/HotspotBalancing.h"

#include "exahype/repositories/Repository.h"
#include "exahype/repositories/RepositoryFactory.h"
#include "exahype/mappings/TimeStepSizeComputation.h"
#include "exahype/mappings/Sending.h"

#include "tarch/Assertions.h"

#include "tarch/logging/CommandLineLogger.h"

#include "tarch/parallel/Node.h"
#include "tarch/parallel/NodePool.h"
#include "tarch/parallel/FCFSNodePoolStrategy.h"

#include "tarch/multicore/Core.h"
#include "tarch/multicore/MulticoreDefinitions.h"

#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/JoinDataBufferPool.h"
#include "peano/parallel/SendReceiveBufferPool.h"
#include "peano/parallel/loadbalancing/Oracle.h"
#include "peano/parallel/loadbalancing/OracleForOnePhaseWithGreedyPartitioning.h"

#include "peano/geometry/Hexahedron.h"

#include "peano/utils/UserInterface.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"

#include "sharedmemoryoracles/OracleForOnePhaseWithGrainSizeSampling.h"
#include "sharedmemoryoracles/OracleForOnePhaseWithShrinkingGrainSize.h"

#ifdef Parallel
#include "mpibalancing/GreedyBalancing.h"
#include "mpibalancing/FairNodePoolStrategy.h"
#endif
#include "exahype/plotters/Plotter.h"

#include "exahype/solvers/LimitingADERDGSolver.h"


tarch::logging::Log exahype::runners::Runner::_log("exahype::runners::Runner");

exahype::runners::Runner::Runner(const Parser& parser) : _parser(parser) {}

exahype::runners::Runner::~Runner() {}

void exahype::runners::Runner::initDistributedMemoryConfiguration() {
  #ifdef Parallel
  std::string configuration = _parser.getMPIConfiguration();
  if (_parser.getMPILoadBalancingType()==Parser::MPILoadBalancingType::Static) {
    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      if (configuration.find( "FCFS" )!=std::string::npos ) {
        tarch::parallel::NodePool::getInstance().setStrategy(
            new tarch::parallel::FCFSNodePoolStrategy()
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on FCFS answering strategy");
      }
      else if (configuration.find( "fair" )!=std::string::npos ) {
        int ranksPerNode = static_cast<int>(exahype::Parser::getValueFromPropertyString(configuration,"ranks_per_node"));
        if (ranksPerNode<=0) {
          logError( "initDistributedMemoryConfiguration()", "please inform fair balancing how many ranks per node you use through value \"ranks_per_node:XXX\". Read value " << ranksPerNode << " is invalid" );
          ranksPerNode = 1;
        }
        if ( ranksPerNode>=tarch::parallel::Node::getInstance().getNumberOfNodes() ) {
          logWarning( "initDistributedMemoryConfiguration()", "value \"ranks_per_node:XXX\" exceeds total rank count. Reset to 1" );
          ranksPerNode = 1;
        }
        tarch::parallel::NodePool::getInstance().setStrategy(
            new mpibalancing::FairNodePoolStrategy(ranksPerNode)
        );
        logInfo("initDistributedMemoryConfiguration()", "load balancing relies on fair answering strategy with " << ranksPerNode << " rank(s) per node") ;
      }
      else {
        logError("initDistributedMemoryConfiguration()", "no valid load balancing answering strategy specified: use FCFS");
        tarch::parallel::NodePool::getInstance().setStrategy(
            new tarch::parallel::FCFSNodePoolStrategy()
        );
      }
    }

    if ( configuration.find( "greedy" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use greedy load balancing without joins (mpibalancing/GreedyBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
        new mpibalancing::GreedyBalancing(
          getCoarsestGridLevelOfAllSolvers(),
          getCoarsestGridLevelOfAllSolvers()+1
        )
      );
    }
    else if ( configuration.find( "hotspot" )!=std::string::npos ) {
      logInfo("initDistributedMemoryConfiguration()", "use global hotspot elimination without joins (mpibalancing/StaticBalancing)");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new mpibalancing::HotspotBalancing(false,getCoarsestGridLevelOfAllSolvers()+1)
      );
    }
    else {
      logError("initDistributedMemoryConfiguration()", "no valid load balancing configured. Use greedy load balancing without joins");
      peano::parallel::loadbalancing::Oracle::getInstance().setOracle(
          new peano::parallel::loadbalancing::OracleForOnePhaseWithGreedyPartitioning(false)
      );
    }
  }
  else {
    logError("initDistributedMemoryConfiguration()", "only MPI static load balancing supported so far. ");
  }

  // end of static load balancing


  tarch::parallel::NodePool::getInstance().restart();
  tarch::parallel::NodePool::getInstance().waitForAllNodesToBecomeIdle();

  tarch::parallel::Node::getInstance().setDeadlockTimeOut(_parser.getMPITimeOut());
  tarch::parallel::Node::getInstance().setTimeOutWarning(_parser.getMPITimeOut()/2);
  logInfo("initDistributedMemoryConfiguration()", "use MPI time out of " << _parser.getMPITimeOut() << " (warn after half the timeout span)");

  const int bufferSize = _parser.getMPIBufferSize();
  peano::parallel::SendReceiveBufferPool::getInstance().setBufferSize(bufferSize);
  peano::parallel::JoinDataBufferPool::getInstance().setBufferSize(bufferSize);
  logInfo("initDistributedMemoryConfiguration()", "use MPI buffer size of " << bufferSize);

  if ( _parser.getSkipReductionInBatchedTimeSteps() ) {
    logInfo("initDistributedMemoryConfiguration()", "allow ranks to skip reduction" );
    exahype::mappings::Sending::SkipReductionInBatchedTimeSteps = true;
  }
  else {
    logWarning("initDistributedMemoryConfiguration()", "ranks are not allowed to skip any reduction (might harm performance). Use optimisation section to switch feature on" );
    exahype::mappings::Sending::SkipReductionInBatchedTimeSteps = false;
  }
  #endif
}


void exahype::runners::Runner::shutdownDistributedMemoryConfiguration() {
#ifdef Parallel
  tarch::parallel::NodePool::getInstance().terminate();
  exahype::repositories::RepositoryFactory::getInstance().shutdownAllParallelDatatypes();
#endif
}

void exahype::runners::Runner::initSharedMemoryConfiguration() {
  #ifdef SharedMemoryParallelisation
  const int numberOfThreads = _parser.getNumberOfThreads();
  tarch::multicore::Core::getInstance().configure(numberOfThreads);

  switch (_parser.getMulticoreOracleType()) {
  case Parser::MulticoreOracleType::Dummy:
    logInfo("initSharedMemoryConfiguration()",
        "use dummy shared memory oracle");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new peano::datatraversal::autotuning::OracleForOnePhaseDummy()
    );
    break;
  case Parser::MulticoreOracleType::Autotuning:
    logInfo("initSharedMemoryConfiguration()",
        "use autotuning shared memory oracle");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithShrinkingGrainSize(
          tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1,
          true
        ));
    break;
  case Parser::MulticoreOracleType::GrainSizeSampling:
    logInfo("initSharedMemoryConfiguration()",
        "use shared memory oracle sampling");
    peano::datatraversal::autotuning::Oracle::getInstance().setOracle(
        new sharedmemoryoracles::OracleForOnePhaseWithGrainSizeSampling(
            64,
            true    // logarithmicDistribution
        ));
    break;
  }

  std::ifstream f(_parser.getMulticorePropertiesFile().c_str());
  if (f.good()) {
    peano::datatraversal::autotuning::Oracle::getInstance().loadStatistics(
        _parser.getMulticorePropertiesFile());
  }
  f.close();
  #endif
}


void exahype::runners::Runner::initDataCompression() {
  exahype::solvers::ADERDGSolver::CompressionAccuracy = _parser.getDoubleCompressionFactor();

  if (exahype::solvers::ADERDGSolver::CompressionAccuracy==0.0) {
    logInfo( "initDataCompression()", "switched off any data compression");
  }
  else {
    if (!_parser.getFuseAlgorithmicSteps()) {
      logError( "initDataCompression()", "data compression is not supported if you don't use the fused time stepping");
      exahype::solvers::ADERDGSolver::CompressionAccuracy = 0.0;
    }
    else {
      exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread = _parser.getSpawnDoubleCompressionAsBackgroundTask();
      logInfo( "initDataCompression()", "store all data with accuracy of " << exahype::solvers::ADERDGSolver::CompressionAccuracy << ". Use background threads for data conversion=" << exahype::solvers::ADERDGSolver::SpawnCompressionAsBackgroundThread);
    }
  }
}


void exahype::runners::Runner::shutdownSharedMemoryConfiguration() {
#ifdef SharedMemoryParallelisation
  switch (_parser.getMulticoreOracleType()) {
  case Parser::MulticoreOracleType::Dummy:
    break;
  case Parser::MulticoreOracleType::Autotuning:
  case Parser::MulticoreOracleType::GrainSizeSampling:
  #ifdef Parallel
    if (tarch::parallel::Node::getInstance().getRank()==tarch::parallel::Node::getInstance().getNumberOfNodes()-1) {
      logInfo("shutdownSharedMemoryConfiguration()",
          "wrote statistics into file " << _parser.getMulticorePropertiesFile()
          << ". Dump from all other ranks subpressed to avoid file races"
      );
      peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
          _parser.getMulticorePropertiesFile());
    }
  #else
    logInfo("shutdownSharedMemoryConfiguration()",
        "wrote statistics into file "
        << _parser.getMulticorePropertiesFile());
    peano::datatraversal::autotuning::Oracle::getInstance().plotStatistics(
        _parser.getMulticorePropertiesFile());
  #endif
    break;
  }
#endif
}


int exahype::runners::Runner::getCoarsestGridLevelOfAllSolvers() const {
  double boundingBox = _parser.getBoundingBoxSize()(0);
  double hMax        = exahype::solvers::Solver::getCoarsestMeshSizeOfAllSolvers();

  int    result      = 1;
  double currenthMax = std::numeric_limits<double>::max();
  while (currenthMax>hMax) {
    currenthMax = boundingBox / threePowI(result);
    result++;
  }

  logDebug( "getCoarsestGridLevelOfAllSolvers()", "regular grid depth of " << result << " creates grid with h_max=" << currenthMax );
  return std::max(3,result);
}


exahype::repositories::Repository* exahype::runners::Runner::createRepository() const {
  // Geometry is static as we need it survive the whole simulation time.
  static peano::geometry::Hexahedron geometry(
      _parser.getDomainSize(),
      _parser.getOffset());

  logDebug(
      "createRepository(...)",
      "create computational domain at " << _parser.getOffset() <<
      " of width/size " << _parser.getDomainSize() <<
      ". bounding box has size " << _parser.getBoundingBoxSize() <<
      ". grid regular up to level " << getCoarsestGridLevelOfAllSolvers() << " (level 1 is coarsest available cell in tree)" );
#ifdef Parallel
  const double boundingBoxScaling = static_cast<double>(getCoarsestGridLevelOfAllSolvers()) / (static_cast<double>(getCoarsestGridLevelOfAllSolvers())-2);
  assertion4(boundingBoxScaling>=1.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize() );
  const double boundingBoxShift   = (1.0-boundingBoxScaling)/2.0;
  assertion5(boundingBoxShift<=0.0, boundingBoxScaling, getCoarsestGridLevelOfAllSolvers(), _parser.getDomainSize(), _parser.getBoundingBoxSize(), boundingBoxScaling );

  logInfo(
      "createRepository(...)",
      "increase domain artificially by " << boundingBoxScaling << " and shift bounding box by " << boundingBoxShift << " to simplify load balancing along boundary");
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _parser.getBoundingBoxSize()*boundingBoxScaling,
      _parser.getOffset()+boundingBoxShift*_parser.getBoundingBoxSize()
  );
#else
  return exahype::repositories::RepositoryFactory::getInstance().createWithSTDStackImplementation(
      geometry,
      _parser.getBoundingBoxSize(),
      _parser.getOffset()
  );
#endif
}


int exahype::runners::Runner::run() {
  exahype::repositories::Repository* repository = createRepository();

  initDistributedMemoryConfiguration();
  initSharedMemoryConfiguration();
  initDataCompression();

  int result = 0;
  if ( _parser.isValid() ) {
    // We have to do this for any rank.
    exahype::State::FuseADERDGPhases = _parser.getFuseAlgorithmicSteps();

    if (tarch::parallel::Node::getInstance().isGlobalMaster()) {
      result = runAsMaster(*repository);
    }
    #ifdef Parallel
    else {
      result = runAsWorker(*repository);
    }
    #endif
  }
  else {
    logError( "run(...)", "do not run code as parser reported errors" );
    result = 1;
  }

  shutdownSharedMemoryConfiguration();
  shutdownDistributedMemoryConfiguration();

  delete repository;

  return result;
}

void exahype::runners::Runner::createGrid(exahype::repositories::Repository& repository) {
#ifdef Parallel
  const bool UseStationaryCriterion = tarch::parallel::Node::getInstance().getNumberOfNodes()==1;
#else
  const bool UseStationaryCriterion = true;
#endif

  int gridSetupIterations = 0;
  repository.switchToMeshRefinement();

  int gridSetupIterationsToRun = 3;
  while (gridSetupIterationsToRun>0) {
    repository.iterate();
    gridSetupIterations++;

    if ( UseStationaryCriterion && repository.getState().isGridStationary() ) {
      gridSetupIterationsToRun--;
    }
    else if ( !repository.getState().isGridBalanced() && tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0 ) {
      gridSetupIterationsToRun=3;  // we need at least 3 sweeps to recover from ongoing balancing
    }
    else if ( !repository.getState().isGridBalanced()  ) {
      gridSetupIterationsToRun=1;  // one additional step to get adjacency right
    }
    else {
      gridSetupIterationsToRun--;
    }

    #if defined(TrackGridStatistics) && defined(Asserts)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #elif defined(Asserts)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", state=" << repository.getState().toString() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #elif defined(TrackGridStatistics)
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", max-level=" << repository.getState().getMaxLevel() <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #else
    logInfo("createGrid()",
        "grid setup iteration #" << gridSetupIterations <<
        ", idle-nodes=" << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
    );
    #endif

    #if !defined(Parallel)
    logInfo("createGrid(...)", "memoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
    #endif

    #ifdef Asserts
    if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
      DataHeap::getInstance().plotStatistics();
      peano::heap::PlainCharHeap::getInstance().plotStatistics();
    }
    #endif
  }

  logInfo("createGrid(Repository)", "finished grid setup after " << gridSetupIterations << " iterations" );

  if (
    tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()>0
    &&
    tarch::parallel::Node::getInstance().getNumberOfNodes()>1
  ) {
    logWarning( "createGrid(Repository)", "there are still " << tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes() << " ranks idle" )
  }

#ifdef Parallel
  // Might be too restrictive for later runs. Remove but keep warning from above
  assertion( tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0 );
#endif
}


int exahype::runners::Runner::runAsMaster(exahype::repositories::Repository& repository) {
  peano::utils::UserInterface::writeHeader();

  if (!exahype::solvers::RegisteredSolvers.empty()) {
    initSolverTimeStepData();
    repository.getState().switchToPreAMRContext();
    createGrid(repository);
    /*
     * Set ADER-DG corrector time stamp and finite volumes time stamp.
     * Compute ADER-DG corrector time step size implicitly and finite volumes time step size.
     * (Implicitly means here that we set the predictor time step size but after the next newTimeStep(...)
     * the corrector time step size is set as this predictor time step size.)
     *
     * Note that it is important that we run SolutionAdjustmentAnd
     * GlobalTimeStepComputation directly after the grid setup
     * since we receive here the metadata
     * that was sent in the last iteration of the grid setup.
     */
    initSolverTimeStepData();

    logInfo( "runAsMaster(...)", "start to initialise all data and to compute first time step size" );

    repository.getState().switchToInitialConditionAndTimeStepSizeComputationContext();
    repository.switchToInitialConditionAndTimeStepSizeComputation();
    repository.iterate();
    logInfo( "runAsMaster(...)", "initialised all data and computed first time step size" );

    // TODO(Dominic):
    if (!exahype::State::fuseADERDGPhases() &&
        exahype::solvers::LimitingADERDGSolver::limiterDomainOfOneSolverHasChanged()) {
      initSolverTimeStepData();

      updateLimiterDomain(repository);
    }

    /*
     * Compute current first predictor based on current time step size.
     * Set current time step size as old time step size of next iteration.
     * Compute the current time step size of the next iteration.
     */
    repository.getState().switchToPredictionAndFusedTimeSteppingInitialisationContext();
    bool plot = exahype::plotters::isAPlotterActive(
        solvers::Solver::getMinSolverTimeStampOfAllSolvers());
    if (plot) {
      #if DIMENSIONS==2
      repository.switchToPredictionAndFusedTimeSteppingInitialisationAndPlot2d();
      #else
      repository.switchToPredictionAndFusedTimeSteppingInitialisationAndPlot();
      #endif
    } else {
      repository.switchToPredictionAndFusedTimeSteppingInitialisation();
    }
    repository.iterate();
    logInfo("runAsMaster(...)","plotted initial solution (if specified) and computed first predictor");

    /*
     * Finally print the initial time step info.
     */
    printTimeStepInfo(-1,repository);

    validateInitialSolverTimeStepData(_parser.getFuseAlgorithmicSteps());

    const double simulationEndTime = _parser.getSimulationEndTime();

    logDebug("runAsMaster(...)","min solver time stamp: "     << solvers::Solver::getMinSolverTimeStampOfAllSolvers());
    logDebug("runAsMaster(...)","min solver time step size: " << solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers());

    while ((solvers::Solver::getMinSolverTimeStampOfAllSolvers() < simulationEndTime) &&
        tarch::la::greater(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
      // TODO(Dominic): This plotting strategy might be an issue if we use LTS.
      // see issue #103
      bool plot = exahype::plotters::isAPlotterActive(
          solvers::Solver::estimateMinNextSolverTimeStampOfAllSolvers());

      if (_parser.getFuseAlgorithmicSteps()) {
        repository.getState().setTimeStepSizeWeightForPredictionRerun(
            _parser.getFuseAlgorithmicStepsFactor());

        int numberOfStepsToRun = 1;
        if (plot) {
          numberOfStepsToRun = 0;
        }
        else if (solvers::Solver::allSolversUseTimeSteppingScheme(solvers::Solver::TimeStepping::GlobalFixed)) {
          /**
           * This computation is optimistic. If we were pessimistic, we had to
           * use the max solver time step size. However, this is not necessary
           * here, as we half the time steps anyway.
           */
          if (solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers()>0.0) {
            const double timeIntervalTillNextPlot = std::min(exahype::plotters::getTimeOfNextPlot(),simulationEndTime) - solvers::Solver::getMaxSolverTimeStampOfAllSolvers();
            numberOfStepsToRun = std::floor( timeIntervalTillNextPlot / solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers() * _parser.getTimestepBatchFactor() );
          }
          numberOfStepsToRun = numberOfStepsToRun<1 ? 1 : numberOfStepsToRun;
        }

        runOneTimeStampWithFusedAlgorithmicSteps(
          repository,
          numberOfStepsToRun,
          _parser.getExchangeBoundaryDataInBatchedTimeSteps() && repository.getState().isGridStationary()
        );
        recomputePredictorIfNecessary(repository);
        printTimeStepInfo(numberOfStepsToRun,repository);
      } else {
        runOneTimeStampWithThreeSeparateAlgorithmicSteps(repository, plot);
      }

      logDebug("runAsMaster(...)", "state=" << repository.getState().toString());
    }
    if ( tarch::la::equals(solvers::Solver::getMinSolverTimeStepSizeOfAllSolvers(), 0.0)) {
      logWarning("runAsMaster(...)","Minimum solver time step size is zero (up to machine precision).");
    }

    repository.logIterationStatistics(false);
  }

  repository.terminate();

  return 0;
}

void exahype::runners::Runner::initSolverTimeStepData() {
  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    p->initSolverTimeStepData(0.0);
  }
}

void exahype::runners::Runner::validateInitialSolverTimeStepData(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0);
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }
      }
      break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        assertionEquals(aderdgSolver->getPreviousMinCorrectorTimeStepSize(),0.0);
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertionEquals(aderdgSolver->getMinCorrectorTimeStamp(),0.0);
        if (fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinPredictorTimeStepSize());
        } else {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),0.0);
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Finite Volumes
        auto* finiteVolumesSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter().get();
        assertionEquals(finiteVolumesSolver->getPreviousMinTimeStepSize(),0.0);
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),0.0);
        assertion1(std::isfinite(finiteVolumesSolver->getMinTimeStepSize()),finiteVolumesSolver->getMinTimeStepSize());
        assertion1(finiteVolumesSolver->getMinTimeStepSize()>0,finiteVolumesSolver->getMinTimeStepSize());
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(finiteVolumesSolver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Compare ADER-DG vs Finite-Volumes
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        assertionEquals(finiteVolumesSolver->getMinTimeStepSize(),aderdgSolver->getMinPredictorTimeStepSize());
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        auto* finiteVolumesSolver = static_cast<exahype::solvers::FiniteVolumesSolver*>(solver);
        assertionEquals(finiteVolumesSolver->getPreviousMinTimeStepSize(),0.0);
        break;
    }
  }
  #endif
}

void exahype::runners::Runner::updateLimiterDomain(exahype::repositories::Repository& repository) {
  logInfo("updateLimiterDomain(...)","start to update limiter domain");
  repository.getState().switchToLimiterStatusSpreadingContext();
  repository.switchToLimiterStatusSpreading();
  repository.iterate();

  /**
   * We need to gather information from all neighbours
   * of a cell before we can determine the unified
   * limiter status value of the cell.
   *
   * We thus need two extra iterations to send and receive
   * the limiter status of remote neighbours.
   */
  #ifdef Parallel
  logInfo("updateLimiterDomain(...)","merge limiter status of remote neighbours");
  repository.switchToLimiterStatusMergingAndSpreadingMPI();
  repository.iterate();
  repository.switchToLimiterStatusMergingMPI();
  repository.iterate();
  #endif

  repository.getState().switchToReinitialisationContext();
  logInfo("updateLimiterDomain(...)","reinitialise cells");
  repository.switchToReinitialisation();
  repository.iterate();

  logInfo("updateLimiterDomain(...)","recompute solution in troubled cells");
  repository.getState().switchToRecomputeSolutionAndTimeStepSizeComputationContext();
  repository.switchToSolutionRecomputationAndTimeStepSizeComputation();
  repository.iterate();

  logInfo("updateLimiterDomain(...)","updated limiter domain");
}

void exahype::runners::Runner::printTimeStepInfo(int numberOfStepsRanSinceLastCall, const exahype::repositories::Repository& repository) {
  double currentMinTimeStamp    = std::numeric_limits<double>::max();
  double currentMinTimeStepSize = std::numeric_limits<double>::max();
  double nextMinTimeStepSize    = std::numeric_limits<double>::max();

  static int n = 0;
  if (numberOfStepsRanSinceLastCall==0) {
    n++;
  }
  else if (numberOfStepsRanSinceLastCall>0) {
    n+=numberOfStepsRanSinceLastCall;
  }

  for (const auto& p : exahype::solvers::RegisteredSolvers) {
    currentMinTimeStamp =
        std::min(currentMinTimeStamp, p->getMinTimeStamp());
    currentMinTimeStepSize =
        std::min(currentMinTimeStepSize, p->getMinTimeStepSize());
    nextMinTimeStepSize =
        std::min(nextMinTimeStepSize, p->getMinNextTimeStepSize());

    #if defined(Debug) || defined(Asserts)
    switch(p->getType()) {
      case exahype::solvers::Solver::Type::ADERDG:
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min         =" << static_cast<exahype::solvers::ADERDGSolver*>(p)->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG correction: dt_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinCorrectorTimeStepSize());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: t_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStamp());
        logInfo("startNewTimeStep(...)",
                "\tADER-DG prediction: dt_min         =" << static_cast<exahype::solvers::LimitingADERDGSolver*>(p)->getSolver()->getMinPredictorTimeStepSize());
        break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
    #endif
  }

  logInfo("startNewTimeStep(...)",
      "step " << n << "\tt_min          =" << currentMinTimeStamp);

  logInfo("startNewTimeStep(...)",
      "\tdt_min         =" << currentMinTimeStepSize);

  #if !defined(Parallel)
  // memory consumption on rank 0 would not make any sense
  logInfo("startNewTimeStep(...)",
      "\tmemoryUsage    =" << peano::utils::UserInterface::getMemoryUsageMB() << " MB");
  #ifdef Asserts
  if (exahype::solvers::ADERDGSolver::CompressionAccuracy>0.0) {
    DataHeap::getInstance().plotStatistics();
    peano::heap::PlainCharHeap::getInstance().plotStatistics();

    logInfo(
      "startNewTimeStep(...)",
      "\tpiped-uncompressed-byes=" << exahype::solvers::ADERDGSolver::PipedUncompressedBytes
      << "\tpiped-compressed-byes=" << exahype::solvers::ADERDGSolver::PipedCompressedBytes
      << "\tcompression-rate=" << (exahype::solvers::ADERDGSolver::PipedCompressedBytes/exahype::solvers::ADERDGSolver::PipedUncompressedBytes)
    );
  }
  #endif

  #if defined(TrackGridStatistics)
  logInfo(
    "startNewTimeStep(...)",
    "\tinner cells/inner unrefined cells=" << repository.getState().getNumberOfInnerCells()
    << "/" << repository.getState().getNumberOfInnerLeafCells() );
  logInfo(
    "startNewTimeStep(...)",
    "\tinner max/min mesh width=" << repository.getState().getMaximumMeshWidth()
    << "/" << repository.getState().getMinimumMeshWidth()
    );
  logInfo(
    "startNewTimeStep(...)",
    "\tmax level=" << repository.getState().getMaxLevel()
    );
  #endif

  #endif
  #if defined(Debug) || defined(Asserts)
  tarch::logging::CommandLineLogger::getInstance().closeOutputStreamAndReopenNewOne();
  #endif
}


void exahype::runners::Runner::runOneTimeStampWithFusedAlgorithmicSteps(
    exahype::repositories::Repository& repository, int numberOfStepsToRun, bool exchangeBoundaryData) {
  /*
   * The adapter below performs the following steps:
   *
   * 1. Exchange the fluctuations using the predictor computed in the previous
   *sweep
   *    and the corrector time stemp size.
   * 2. Perform the corrector step using the corrector update and the corrector
   *time step size.
   *    This is a cell-local operation. Thus we immediately obtain the
   *cell-local current solution.
   * 3. Perform the predictor step using the cell-local current solution and the
   *predictor time step size.
   * 4. Compute the cell-local time step sizes
   */
  repository.getState().switchToADERDGTimeStepContext();
  if (numberOfStepsToRun==0) {
    repository.switchToPlotAndADERDGTimeStep();
    repository.iterate();
  } else {
    repository.switchToADERDGTimeStep();
    repository.iterate(numberOfStepsToRun,exchangeBoundaryData);
  }

  // reduction/broadcast barrier
}

void exahype::runners::Runner::recomputePredictorIfNecessary(
    exahype::repositories::Repository& repository) {
  // Must be evaluated before we start a new time step
  //  bool stabilityConditionWasHarmed = setStableTimeStepSizesIfStabilityConditionWasHarmed(factor);
  // Note that it is important to switch the time step sizes, i.e,
  // start a new time step, before we recompute the predictor.

  if (repository.getState().stabilityConditionOfOneSolverWasViolated()) {
    logInfo("startNewTimeStep(...)",
        "\t\t Space-time predictor must be recomputed.");

    repository.getState().switchToPredictionRerunContext();
    repository.switchToPredictionRerun();
    repository.iterate();
  }
}

void exahype::runners::Runner::runOneTimeStampWithThreeSeparateAlgorithmicSteps(
    exahype::repositories::Repository& repository, bool plot) {
  // Only one time step (predictor vs. corrector) is used in this case.
  logInfo("runOneTimeStampWithThreeSeparateAlgorithmicSteps(...)","merge neighbours");

  repository.getState().switchToNeighbourDataMergingContext();
  repository.switchToNeighbourDataMerging();  // Riemann -> face2face
  repository.iterate(); // todo uncomment

  logInfo("runOneTimeStampWithThreeSeparateAlgorithmicSteps(...)","update solution");

  repository.getState().switchToSolutionUpdateContext();
  repository.switchToSolutionUpdate();  // Face to cell + Inside cell
  repository.iterate();

  logInfo("runOneTimeStampWithThreeSeparateAlgorithmicSteps(...)","compute new time step size");

  repository.getState().switchToTimeStepSizeComputationContext();
  repository.switchToTimeStepSizeComputation();
  repository.iterate();

  // We mimic the flow of the fused time stepping scheme here.
  // Updating the limiter domain is thus done after the time step
  // size computation.
  if (exahype::solvers::LimitingADERDGSolver::limiterDomainOfOneSolverHasChanged()) {
    updateLimiterDomain(repository);
  }

  printTimeStepInfo(1,repository);

  /*
   * Compute current first predictor based on current time step size.
   * Set current time step size as old time step size of next iteration.
   * Compute the current time step size of the next iteration.
   *
   * TODO(Dominic): Limiting: There is an issue with the prediction in
   * the limiting context. Since we overwrite the update here again.
   * A rollback is thus not possible anymore.
   * The only way out of here would be to store an old and new
   * ADER-DG solution similar to the finite volumes solver. This is the reason
   * why we currently only offer the limiting for
   * the non-fused time stepping variant.
   */
  repository.getState().switchToPredictionContext();
  if (plot) {
    #if DIMENSIONS==2
    repository.switchToPredictionAndPlot2d();
    #else
    repository.switchToPredictionAndPlot();
    #endif
  } else {
    repository.switchToPrediction();   // Cell onto faces
  }
  repository.iterate();
}

void exahype::runners::Runner::validateSolverTimeStepDataForThreeAlgorithmicPhases(const bool fuseADERDGPhases) const {
  #ifdef Asserts
  for (auto* solver : exahype::solvers::RegisteredSolvers) {
    assertionEquals(solver->getMinTimeStamp(),0.0);
    assertion1(std::isfinite(solver->getMinTimeStepSize()),solver->getMinTimeStepSize());
    assertion1(solver->getMinTimeStepSize()>0,solver->getMinTimeStepSize());

    switch(solver->getTimeStepping()) {
      case exahype::solvers::Solver::TimeStepping::Global:
        assertionEquals(solver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
        break;
      case exahype::solvers::Solver::TimeStepping::GlobalFixed:
        break;
    }
    switch (solver->getType()) {
      case exahype::solvers::Solver::Type::ADERDG: {
        auto* aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

      }  break;
      case exahype::solvers::Solver::Type::LimitingADERDG: {
        // ADER-DG
        auto* aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        assertion1(std::isfinite(aderdgSolver->getMinPredictorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(std::isfinite(aderdgSolver->getMinCorrectorTimeStepSize()),aderdgSolver->getMinPredictorTimeStepSize());
        assertion1(aderdgSolver->getMinCorrectorTimeStamp() > 0.0,aderdgSolver->getMinCorrectorTimeStamp());
        if (!fuseADERDGPhases) {
          assertionEquals(aderdgSolver->getMinPredictorTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        }
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(aderdgSolver->getMinNextPredictorTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Finite Volumes
        auto* finiteVolumesSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getLimiter().get();
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),0.0);
        assertion1(std::isfinite(finiteVolumesSolver->getMinTimeStepSize()),finiteVolumesSolver->getMinTimeStepSize());
        assertion1(finiteVolumesSolver->getMinTimeStepSize()>0,finiteVolumesSolver->getMinTimeStepSize());
        switch(solver->getTimeStepping()) {
          case exahype::solvers::Solver::TimeStepping::Global:
            assertionEquals(finiteVolumesSolver->getMinNextTimeStepSize(),std::numeric_limits<double>::max());
            break;
          case exahype::solvers::Solver::TimeStepping::GlobalFixed:
            break;
        }

        // Compare ADER-DG vs Finite-Volumes
        assertionEquals(finiteVolumesSolver->getMinTimeStamp(),aderdgSolver->getMinCorrectorTimeStamp());
        assertionEquals(finiteVolumesSolver->getMinTimeStepSize(),aderdgSolver->getMinPredictorTimeStepSize());
      } break;
      case exahype::solvers::Solver::Type::FiniteVolumes:
        break;
    }
  }
  #endif
}
