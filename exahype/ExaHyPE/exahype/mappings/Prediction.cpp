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
 
#include "exahype/mappings/Prediction.h"

#include "peano/datatraversal/autotuning/Oracle.h"
#include "peano/utils/Globals.h"

#include "tarch/multicore/Loop.h"
#include "tarch/multicore/Lock.h"

#include "peano/utils/Loop.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "exahype/solvers/ADERDGSolver.h"
#include "exahype/solvers/LimitingADERDGSolver.h"

#include "exahype/amr/AdaptiveMeshRefinement.h"

#include "peano/utils/UserInterface.h"

#include <algorithm>
#include <mm_malloc.h> //g++
#include <cstring> //memset

peano::CommunicationSpecification
exahype::mappings::Prediction::communicationSpecification() {
  return peano::CommunicationSpecification(
      peano::CommunicationSpecification::ExchangeMasterWorkerData::MaskOutMasterWorkerDataAndStateExchange,
      peano::CommunicationSpecification::ExchangeWorkerMasterData::MaskOutWorkerMasterDataAndStateExchange,
      true);
}

peano::MappingSpecification
exahype::mappings::Prediction::enterCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::WholeTree,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Prediction::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

// The remainder specs all are nop
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Prediction::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}
peano::MappingSpecification
exahype::mappings::Prediction::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}
peano::MappingSpecification
exahype::mappings::Prediction::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidCoarseGridRaces,true);
}

tarch::logging::Log exahype::mappings::Prediction::_log(
    "exahype::mappings::Prediction");


void exahype::mappings::Prediction::prepareTemporaryVariables() {
  assertion(_tempSpaceTimeUnknowns    ==nullptr);
  assertion(_tempSpaceTimeFluxUnknowns==nullptr);
  assertion(_tempUnknowns             ==nullptr);
  assertion(_tempFluxUnknowns         ==nullptr);
  assertion(_tempStateSizedVectors    ==nullptr);
  assertion(_tempPointForceSources    ==nullptr);

  int numberOfSolvers        = exahype::solvers::RegisteredSolvers.size();
  _tempSpaceTimeUnknowns     = new double**[numberOfSolvers]; // == lQi, lQi_old, rhs, rhs_0 (unchanged by optimisation)
  _tempSpaceTimeFluxUnknowns = new double**[numberOfSolvers]; // == lFi, gradQ
  _tempUnknowns              = new double* [numberOfSolvers]; // == lQhi
  _tempFluxUnknowns          = new double* [numberOfSolvers]; // == lFhi
  _tempStateSizedVectors     = new double* [numberOfSolvers]; // == BGradQ
  _tempPointForceSources     = new double* [numberOfSolvers];

  exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;

  int solverNumber=0;
  for (auto solver : exahype::solvers::RegisteredSolvers) {
    switch( solver->getType() ) {
    case exahype::solvers::Solver::Type::ADERDG:
      aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
      break;
    case exahype::solvers::Solver::Type::LimitingADERDG:
      aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
      break;
    default:
      aderdgSolver = nullptr;
      break;
    }

    if (aderdgSolver!=nullptr) {
      if(aderdgSolver->alignTempArray()) {
        _tempSpaceTimeUnknowns[solverNumber] = new double*[4];
        for (int i=0; i<4; ++i) { // max; see spaceTimePredictorNonlinear
          _tempSpaceTimeUnknowns[solverNumber][i] =
              (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize(), ALIGNMENT);
          std::memset(_tempSpaceTimeUnknowns[solverNumber][i], 0, sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize());
        }
        //
        _tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          _tempSpaceTimeFluxUnknowns[solverNumber][i] =
              (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeFluxUnknownsSize(), ALIGNMENT);
          std::memset(_tempSpaceTimeFluxUnknowns[solverNumber][i], 0, sizeof(double)*aderdgSolver->getTempSpaceTimeFluxUnknownsSize());
        }
        //
        _tempUnknowns    [solverNumber]      = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempUnknownsSize(), ALIGNMENT);
        //
        _tempFluxUnknowns[solverNumber]      = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempFluxUnknownsSize(), ALIGNMENT);
         //
        _tempStateSizedVectors[solverNumber] = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getNumberOfVariables(), ALIGNMENT);
        
        if(aderdgSolver->isDummyKRequired()) { //TODO KD
           _tempPointForceSources    [solverNumber] = (double *) _mm_malloc(sizeof(double)*aderdgSolver->getTempSpaceTimeUnknownsSize(), ALIGNMENT);
        } else {
           _tempPointForceSources    [solverNumber] = nullptr;
        }
      } else {
        _tempSpaceTimeUnknowns[solverNumber] = new double*[4];
        for (int i=0; i<4; ++i) { // max; see spaceTimePredictorNonlinear
          _tempSpaceTimeUnknowns[solverNumber][i] =
              new double[aderdgSolver->getTempSpaceTimeUnknownsSize()]();
        }
        //
        _tempSpaceTimeFluxUnknowns[solverNumber] = new double*[2];
        for (int i=0; i<2; ++i) { // max; see spaceTimePredictorNonlinear
          _tempSpaceTimeFluxUnknowns[solverNumber][i] =
              new double[aderdgSolver->getTempSpaceTimeFluxUnknownsSize()]();
        }
        //
        _tempUnknowns    [solverNumber]      = new double[aderdgSolver->getTempUnknownsSize()]; 
        //
        _tempFluxUnknowns[solverNumber]      = new double[aderdgSolver->getTempFluxUnknownsSize()]; 
         //
        _tempStateSizedVectors[solverNumber] = new double[aderdgSolver->getNumberOfVariables()]; 
        
        if(aderdgSolver->isDummyKRequired()) { //TODO KD
           _tempPointForceSources    [solverNumber] = new double[aderdgSolver->getTempSpaceTimeUnknownsSize()];
        } else {
           _tempPointForceSources    [solverNumber] = nullptr;
        } 
      }      
    } else {
      _tempSpaceTimeUnknowns    [solverNumber] = nullptr;
      _tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
      _tempUnknowns             [solverNumber] = nullptr;
      _tempFluxUnknowns         [solverNumber] = nullptr;
      _tempStateSizedVectors    [solverNumber] = nullptr;
      _tempPointForceSources    [solverNumber] = nullptr;
    }

    ++solverNumber;
  }
}

void exahype::mappings::Prediction::deleteTemporaryVariables() {
  if (_tempSpaceTimeUnknowns!=nullptr) {
    assertion(_tempSpaceTimeFluxUnknowns!=nullptr);
    assertion(_tempUnknowns             !=nullptr);
    assertion(_tempFluxUnknowns         !=nullptr);
    assertion(_tempStateSizedVectors    !=nullptr);
    assertion(_tempPointForceSources    !=nullptr);

    int solverNumber=0;
    exahype::solvers::ADERDGSolver* aderdgSolver = nullptr;
    
    for (auto solver : exahype::solvers::RegisteredSolvers) {
      switch( solver->getType() ) {
      case exahype::solvers::Solver::Type::ADERDG:
        aderdgSolver = static_cast<exahype::solvers::ADERDGSolver*>(solver);
        break;
      case exahype::solvers::Solver::Type::LimitingADERDG:
        aderdgSolver = static_cast<exahype::solvers::LimitingADERDGSolver*>(solver)->getSolver().get();
        break;
      default:
        aderdgSolver = nullptr;
        break;
      }
      
      if (aderdgSolver!=nullptr) {
        if(aderdgSolver->alignTempArray()) {
          //
          for (int i=0; i<4; ++i) {
            _mm_free(_tempSpaceTimeUnknowns[solverNumber][i]);
          }
          delete[] _tempSpaceTimeUnknowns[solverNumber];
          _tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            _mm_free(_tempSpaceTimeFluxUnknowns[solverNumber][i]);
          }
          delete[] _tempSpaceTimeFluxUnknowns[solverNumber];
          _tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          _mm_free(_tempUnknowns[solverNumber]);
          _tempUnknowns[solverNumber] = nullptr;
          //
          _mm_free(_tempFluxUnknowns[solverNumber]);
          _tempFluxUnknowns[solverNumber] = nullptr;
          //
          _mm_free(_tempStateSizedVectors[solverNumber]);
          _tempStateSizedVectors[solverNumber] = nullptr;
          
          if(aderdgSolver->isDummyKRequired()) { //TODO KD
            _mm_free(_tempPointForceSources[solverNumber]);
            _tempPointForceSources[solverNumber] = nullptr;
          }
        } else {
          //
          for (int i=0; i<4; ++i) {
            delete[] _tempSpaceTimeUnknowns[solverNumber][i];
          }
          delete[] _tempSpaceTimeUnknowns[solverNumber];
          _tempSpaceTimeUnknowns[solverNumber] = nullptr;
          //
          for (int i=0; i<2; ++i) {
            delete[] _tempSpaceTimeFluxUnknowns[solverNumber][i];
          }
          delete[] _tempSpaceTimeFluxUnknowns[solverNumber];
          _tempSpaceTimeFluxUnknowns[solverNumber] = nullptr;
          //
          delete[] _tempUnknowns[solverNumber];
          _tempUnknowns[solverNumber] = nullptr;
          //
          delete[] _tempFluxUnknowns[solverNumber];
          _tempFluxUnknowns[solverNumber] = nullptr;
          //
          delete[] _tempStateSizedVectors[solverNumber];
          _tempStateSizedVectors[solverNumber] = nullptr;
          
          if(aderdgSolver->isDummyKRequired()) { //TODO KD
            delete[] _tempPointForceSources[solverNumber];
            _tempPointForceSources[solverNumber] = nullptr;
          }
        }
        
      }

      ++solverNumber;
    }

    delete[] _tempSpaceTimeUnknowns;
    delete[] _tempSpaceTimeFluxUnknowns;
    delete[] _tempUnknowns;
    delete[] _tempFluxUnknowns;
    delete[] _tempStateSizedVectors;
    delete[] _tempPointForceSources;
    _tempSpaceTimeUnknowns     = nullptr;
    _tempSpaceTimeFluxUnknowns = nullptr;
    _tempUnknowns              = nullptr;
    _tempFluxUnknowns          = nullptr;
    _tempStateSizedVectors     = nullptr;
    _tempPointForceSources     = nullptr;
  }
}

exahype::mappings::Prediction::Prediction() :
        _tempSpaceTimeUnknowns(nullptr),
        _tempSpaceTimeFluxUnknowns(nullptr),
        _tempUnknowns(nullptr),
        _tempFluxUnknowns(nullptr),
        _tempStateSizedVectors(nullptr),
        _tempPointForceSources(nullptr){}

exahype::mappings::Prediction::~Prediction() {
  deleteTemporaryVariables();
}

#if defined(SharedMemoryParallelisation)
exahype::mappings::Prediction::Prediction(
    const Prediction& masterThread) :
  _tempSpaceTimeUnknowns(nullptr),
  _tempSpaceTimeFluxUnknowns(nullptr),
  _tempUnknowns(nullptr),
  _tempFluxUnknowns(nullptr),
  _tempStateSizedVectors(nullptr),
  _tempPointForceSources(nullptr)  {
  prepareTemporaryVariables();
}

void exahype::mappings::Prediction::mergeWithWorkerThread(
    const Prediction& workerThread) {
}
#endif

void exahype::mappings::Prediction::beginIteration(
    exahype::State& solverState) {
  prepareTemporaryVariables();
}

void exahype::mappings::Prediction::endIteration(
    exahype::State& solverState) {
  deleteTemporaryVariables();
}

void exahype::mappings::Prediction::performPredictionAndVolumeIntegral(
                                        exahype::solvers::ADERDGSolver* solver,
                                        exahype::solvers::ADERDGSolver::CellDescription& cellDescription,
                                        exahype::Vertex* const fineGridVertices,
                                        const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) {
  if (cellDescription.getType()==exahype::records::ADERDGCellDescription::Cell) {
    assertion1(cellDescription.getRefinementEvent()==exahype::records::ADERDGCellDescription::None,cellDescription.toString());

//    solver->validateNoNansInADERDGSolver(cellDescription,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[pre]");

    solver->performPredictionAndVolumeIntegral(
        cellDescription,
        _tempSpaceTimeUnknowns    [cellDescription.getSolverNumber()],
        _tempSpaceTimeFluxUnknowns[cellDescription.getSolverNumber()],
        _tempUnknowns             [cellDescription.getSolverNumber()],
        _tempFluxUnknowns         [cellDescription.getSolverNumber()],
        _tempStateSizedVectors    [cellDescription.getSolverNumber()],
        _tempPointForceSources    [cellDescription.getSolverNumber()]);

    solver->validateNoNansInADERDGSolver(cellDescription,fineGridVerticesEnumerator,"exahype::mappings::Prediction::enterCell[post]");
  }
}

void exahype::mappings::Prediction::enterCell(
    exahype::Cell& fineGridCell,
    exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  logTraceInWith4Arguments("enterCell(...)", fineGridCell,
                           fineGridVerticesEnumerator.toString(),
                           coarseGridCell, fineGridPositionOfCell);

  if (fineGridCell.isInitialised()) {
    const int numberOfADERDGCellDescriptions = static_cast<int>(
        exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
            fineGridCell.getCellDescriptionsIndex()).size());
    if (numberOfADERDGCellDescriptions>0) {
      auto grainSize = peano::datatraversal::autotuning::Oracle::getInstance().parallelise(
          numberOfADERDGCellDescriptions, peano::datatraversal::autotuning::MethodTrace::UserDefined9);
      pfor(i, 0, numberOfADERDGCellDescriptions, grainSize.getGrainSize())
        auto& cellDescription = exahype::solvers::ADERDGSolver::getCellDescription(
            fineGridCell.getCellDescriptionsIndex(),i);

        switch (exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]->getType()) {
          case exahype::solvers::Solver::Type::ADERDG: {
            exahype::solvers::ADERDGSolver* solver = static_cast<exahype::solvers::ADERDGSolver*>(
                exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
            solver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),i); // Time step synchr. might be done multiple times per traversal; but this is no issue.
            solver->prepareNextNeighbourMerging(
                fineGridCell.getCellDescriptionsIndex(),i,
                fineGridVertices,fineGridVerticesEnumerator);

            performPredictionAndVolumeIntegral(solver,cellDescription,fineGridVertices,fineGridVerticesEnumerator);
          } break;
          case exahype::solvers::Solver::Type::LimitingADERDG: {
            exahype::solvers::LimitingADERDGSolver* solver = static_cast<exahype::solvers::LimitingADERDGSolver*>(
                exahype::solvers::RegisteredSolvers[cellDescription.getSolverNumber()]);
            solver->synchroniseTimeStepping(fineGridCell.getCellDescriptionsIndex(),i); // Time step synchr. might be done multiple times per traversal; but this is no issue.
            solver->prepareNextNeighbourMerging(
                fineGridCell.getCellDescriptionsIndex(),i,
                fineGridVertices,fineGridVerticesEnumerator);

            if (cellDescription.getLimiterStatus()==exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::Ok
                || cellDescription.getLimiterStatus()==exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsTroubledCell
                || cellDescription.getLimiterStatus()==exahype::solvers::ADERDGSolver::CellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell) {
              performPredictionAndVolumeIntegral(solver->getSolver().get(),cellDescription,fineGridVertices,fineGridVerticesEnumerator);
            }
          } break;
          default:
            break;
        }
      endpfor
      grainSize.parallelSectionHasTerminated();
    }
  }
  logTraceOutWith1Argument("enterCell(...)", fineGridCell);
}
  // TODO(Dominic): Add getters


//
// Below all methods are nop.
//
// ====================================

#ifdef Parallel
void exahype::mappings::Prediction::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

bool exahype::mappings::Prediction::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  // do nothing
  return false;
}

void exahype::mappings::Prediction::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Vertex& localVertex, const exahype::Vertex& masterOrWorkerVertex,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithRemoteDataDueToForkOrJoin(
    exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
    int fromRank, const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithMaster(
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
    exahype::State& masterState) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {
  // do nothing
}

void exahype::mappings::Prediction::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
  // do nothing
}
#endif

void exahype::mappings::Prediction::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}


void exahype::mappings::Prediction::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
  // do nothing
}

void exahype::mappings::Prediction::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
  // do nothing
}

void exahype::mappings::Prediction::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}

void exahype::mappings::Prediction::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {
  // do nothing
}
