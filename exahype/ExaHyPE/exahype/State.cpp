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
 
#include "exahype/State.h"
#include "exahype/Cell.h"
#include "exahype/Vertex.h"

#include "peano/grid/Checkpoint.h"

#include "exahype/solvers/Solver.h"

#include <limits>

bool exahype::State::FuseADERDGPhases = false;

exahype::State::State() : Base() {
  #ifdef Parallel
  _stateData.setGridConstructionState( exahype::records::State::Default );
  _idleRanksAtLastLookup = 0;
  #endif
}

exahype::State::State(const Base::PersistentState& argument) : Base(argument) {
  // do nothing
}

void exahype::State::merge(const exahype::State& anotherState) {
  // do nothing
}

void exahype::State::writeToCheckpoint(
    peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) const {
  // do nothing
}

void exahype::State::readFromCheckpoint(
    const peano::grid::Checkpoint<exahype::Vertex, exahype::Cell>& checkpoint) {
  // do nothing
}

void exahype::State::updateRegularInitialGridRefinementStrategy() {
  assertion( tarch::parallel::Node::getInstance().isGlobalMaster() );

  #ifdef Parallel
  if (
    tarch::parallel::Node::getInstance().getNumberOfNodes()==1
    ||
    _stateData.getGridConstructionState()==exahype::records::State::Aggressive
  ) {
    _stateData.setGridConstructionState( exahype::records::State::Aggressive );
  }
  else if (
    tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()==0
    &&
    isGridStationary()
  ) {
    _stateData.setGridConstructionState( exahype::records::State::Aggressive );
  }
  else if (
    isInvolvedInJoinOrFork()
    ||
    !isTraversalInverted()
    ||
    !isGridStationary()
  ) {
    _stateData.setGridConstructionState( exahype::records::State::Veto );
  }
  else if (
    _idleRanksAtLastLookup != tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes()
  ) {
    _stateData.setGridConstructionState( exahype::records::State::Veto );
    _idleRanksAtLastLookup = tarch::parallel::NodePool::getInstance().getNumberOfIdleNodes();
  }
  else {
    _stateData.setGridConstructionState( exahype::records::State::Default );
  }

  #endif
}

bool exahype::State::refineInitialGridInCreationalEvents() const {
  #ifdef Parallel
  assertion1( !isInvolvedInJoinOrFork() || _stateData.getGridConstructionState() != exahype::records::State::Aggressive, toString() );
  return _stateData.getGridConstructionState() == exahype::records::State::Aggressive;
  #else
  return true;
  #endif
}

bool exahype::State::refineInitialGridInTouchVertexLastTime() const {
  #ifdef Parallel
  return _stateData.getGridConstructionState() != exahype::records::State::Veto;
  #else
  return false;
  #endif
}
