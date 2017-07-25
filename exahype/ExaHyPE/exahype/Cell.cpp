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
 
#include "exahype/Cell.h"
#include "exahype/State.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

#include "tarch/la/ScalarOperations.h"

#include "peano/utils/Loop.h"

#include "kernels/KernelCalls.h"

#include "exahype/records/ADERDGCellDescription.h"

tarch::logging::Log exahype::Cell::_log("exahype::Cell");

exahype::Cell::Cell() : Base() {
  // We initialise cells which are not touched by the
  // createCell(...) events of Peano's spacetree traversal automaton
  // with default ("do-nothing") values.
  _cellData.setCellDescriptionsIndex(
      multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

exahype::Cell::Cell(const Base::DoNotCallStandardConstructor& value)
    : Base(value) {
}

exahype::Cell::Cell(const Base::PersistentCell& argument) : Base(argument) {
  // This constructor is used to create a cell from persistent data.
  // Do not use it. This would overwrite persistent data.
}

bool exahype::Cell::isFaceInside(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  dfor2(v) // Loop over vertices.
    if (v(d) == f && verticesAroundCell[ verticesEnumerator(v) ].isInside()) {
      return true;
    }
  enddforx // v
  return false;
}

#ifdef Parallel
bool exahype::Cell::isAdjacentToRemoteRank(
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  tarch::la::Vector<DIMENSIONS,int> center(1);
  dfor2(v) // Loop over vertices.
    if (verticesAroundCell[ verticesEnumerator(v) ].isAdjacentToRemoteRank()) {
      dfor2(a) // Loop over adjacent ranks. Does also include own rank.
        if (tarch::la::countEqualEntries(v+a,center)==DIMENSIONS-1 && // offset in one direction from center=>face neighbour
            verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
            tarch::parallel::Node::getInstance().getRank()) {
          return true;
        }
      enddforx //a
    }
  enddforx // v

  return false;
}

int exahype::Cell::countListingsOfRemoteRankByInsideVerticesAtFace(
    const int faceIndex,
    exahype::Vertex* const verticesAroundCell,
    const peano::grid::VertexEnumerator& verticesEnumerator) {
  int result = 0;

  const int f = faceIndex % 2;   // "0" indicates a left face, "1" indicates a right face.
  const int d = (faceIndex-f)/2; // The normal direction: 0: x, 1: y, 1: z.

  tarch::la::Vector<DIMENSIONS,int> pos(1); // This is now the center, i.e., (1,1,...,1).
  pos(d) = 2*f;                             // This is a shift from the center by one unit in direction d.

  int faceNeighbourRank = -1; // This variable is introduced to make sure that the adjacent remote rank is unique.
  // TODO(Dominic): Uniqueness is probably guaranteed by the SFC based DD.
  dfor2(v) // Loop over vertices.
    if (verticesAroundCell[ verticesEnumerator(v) ].isAdjacentToRemoteRank()) {
      dfor2(a) // Loop over adjacent ranks. Does also include own rank.
        if (tarch::la::equals(v+a,pos) &&
            verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]!=
            tarch::parallel::Node::getInstance().getRank()) {
          // Increment
          if (faceNeighbourRank==-1) {
            faceNeighbourRank = verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar];
          }
          if (verticesAroundCell[ verticesEnumerator(v) ].getAdjacentRanks()[aScalar]==faceNeighbourRank) {
            result++;
          }
        }
      enddforx // a
    }
  enddforx // v

  // result must be either no connection, edge connection, or whole face connection.
  assertion2(result==0||result==TWO_POWER_D_DIVIDED_BY_TWO/2||result==TWO_POWER_D_DIVIDED_BY_TWO,result,faceIndex);

  return result;
}
#endif

void exahype::Cell::setupMetaData() {
  assertion1(!exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),toString());

  const int cellDescriptionIndex = exahype::solvers::ADERDGSolver::Heap::getInstance().createData(0, 0);
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().createDataForIndex(cellDescriptionIndex,0,0);

  _cellData.setCellDescriptionsIndex(cellDescriptionIndex);
}

void exahype::Cell::shutdownMetaData() {
  assertion1(
    exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  exahype::solvers::ADERDGSolver::Heap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());
  exahype::solvers::FiniteVolumesSolver::Heap::getInstance().deleteData(_cellData.getCellDescriptionsIndex());

  _cellData.setCellDescriptionsIndex(multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex);
}

bool exahype::Cell::isInitialised() const {
  if (_cellData.getCellDescriptionsIndex()!=multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    assertion1( exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
            _cellData.getCellDescriptionsIndex());
    assertion1( exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(_cellData.getCellDescriptionsIndex()),
        _cellData.getCellDescriptionsIndex());
  }  // Dead code elimination will get rid of this loop if Asserts flag is not set.

  return _cellData.getCellDescriptionsIndex()!=multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex;
}

int exahype::Cell::getCellDescriptionsIndex() const {
  return _cellData.getCellDescriptionsIndex();
}

void exahype::Cell::setCellDescriptionsIndex(int cellDescriptionsIndex) {
  _cellData.setCellDescriptionsIndex(cellDescriptionsIndex);
}

void exahype::Cell::addNewCellDescription(
    const int solverNumber,
    const exahype::records::FiniteVolumesCellDescription::Type cellType,
/*
    const exahype::records::FiniteVolumesCellDescription::RefinementEvent
        refinementEvent,
*/
    const int level,
    const int parentIndex,
    const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
    const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    exahype::solvers::FiniteVolumesSolver::Heap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

//  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
//  assertion(solver->getType()==exahype::solvers::Solver::Type::FiniteVolumes);

  exahype::records::FiniteVolumesCellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setLevel(level);
  //newCellDescription.setRefinementEvent(refinementEvent);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Initialise helper variables
//  newCellDescription.setHelperCellNeedsToStoreFaceData(false); // TODO(Dominic): Add to FV cell descr.

  // Default field data indices
  newCellDescription.setSolution(-1);
  newCellDescription.setPreviousSolution(-1);

  exahype::solvers::FiniteVolumesSolver::Heap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);

}


void exahype::Cell::addNewCellDescription(
  const int                                     solverNumber,
  const exahype::records::ADERDGCellDescription::Type cellType,
  const exahype::records::ADERDGCellDescription::RefinementEvent refinementEvent,
  const int                                     level,
  const int                                     parentIndex,
  const tarch::la::Vector<DIMENSIONS, double>&  cellSize,
  const tarch::la::Vector<DIMENSIONS, double>&  cellOffset) {
  if (_cellData.getCellDescriptionsIndex() == multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex) {
    setupMetaData();
  }

  assertion1(
    exahype::solvers::ADERDGSolver::Heap::getInstance().isValidIndex(
      _cellData.getCellDescriptionsIndex()
    ),
    toString());

  assertion2(parentIndex == -1 ||
             parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(parentIndex != _cellData.getCellDescriptionsIndex(),
             parentIndex, _cellData.getCellDescriptionsIndex());

  assertion2(static_cast<unsigned int>(solverNumber) <
                 solvers::RegisteredSolvers.size(),
             solverNumber, exahype::solvers::RegisteredSolvers.size());

//  const solvers::Solver* solver = solvers::RegisteredSolvers[solverNumber];
//  assertion(solver->getType()==exahype::solvers::Solver::Type::ADER_DG);

  exahype::records::ADERDGCellDescription newCellDescription;
  newCellDescription.setSolverNumber(solverNumber);

  // Default AMR settings
  newCellDescription.setType(cellType);
  newCellDescription.setParentIndex(parentIndex);
  newCellDescription.setLevel(level);
  newCellDescription.setRefinementEvent(refinementEvent);

  std::bitset<DIMENSIONS_TIMES_TWO>
      riemannSolvePerformed;  // default construction: no bit set
  newCellDescription.setRiemannSolvePerformed(riemannSolvePerformed);

  // Pass geometry information to the cellDescription description
  newCellDescription.setSize(cellSize);
  newCellDescription.setOffset(cellOffset);

  // Initialise MPI helper variables
  #ifdef Parallel
  newCellDescription.setHasToHoldDataForNeighbourCommunication(false);
  newCellDescription.setHasToHoldDataForMasterWorkerCommunication(false);
  for (int faceIndex = 0; faceIndex < DIMENSIONS_TIMES_TWO; faceIndex++) {
    newCellDescription.setFaceDataExchangeCounter(faceIndex,TWO_POWER_D);
  }
  #endif

  // Default field data indices
  newCellDescription.setSolution(-1);
  newCellDescription.setUpdate(-1);
  newCellDescription.setExtrapolatedPredictor(-1);
  newCellDescription.setFluctuation(-1);

  // Limiter meta data (oscillations identificator)
  newCellDescription.setLimiterStatus(exahype::records::ADERDGCellDescription::LimiterStatus::Ok); // implicit conversion.
  newCellDescription.setSolutionMin(-1);
  newCellDescription.setSolutionMax(-1);

  // Compression
  newCellDescription.setCompressionState(exahype::records::ADERDGCellDescription::CompressionState::Uncompressed);
  newCellDescription.setSolutionAverages(-1);
  newCellDescription.setUpdateAverages(-1);
  newCellDescription.setExtrapolatedPredictorAverages(-1);
  newCellDescription.setFluctuationAverages(-1);

  newCellDescription.setSolutionCompressed(-1);
  newCellDescription.setUpdateCompressed(-1);
  newCellDescription.setExtrapolatedPredictorCompressed(-1);
  newCellDescription.setFluctuationCompressed(-1);

  exahype::solvers::ADERDGSolver::Heap::getInstance()
      .getData(_cellData.getCellDescriptionsIndex())
      .push_back(newCellDescription);
}

int exahype::Cell::getNumberOfADERDGCellDescriptions() const {
  return exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}


int exahype::Cell::getNumberOfFiniteVolumeCellDescriptions() const {
  return exahype::solvers::FiniteVolumesSolver::Heap::getInstance().getData(
      getCellDescriptionsIndex()).size();
}
