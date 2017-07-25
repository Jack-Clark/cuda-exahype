#include "peano/datatraversal/autotuning/OracleForOnePhaseDummy.h"
#include "peano/datatraversal/autotuning/MethodTrace.h"
#include "tarch/Assertions.h"

#include <cstdlib>
#include <limits>
#include <iostream>
#include <fstream>


tarch::logging::Log  peano::datatraversal::autotuning::OracleForOnePhaseDummy::_log( "peano::datatraversal::autotuning::OracleForOnePhaseDummy" );


peano::datatraversal::autotuning::OracleForOnePhaseDummy::OracleForOnePhaseDummy(
  bool useMultithreading                  ,
  int  grainSizeOfUserDefinedRegions      ,
  SplitVertexReadsOnRegularSubtree  splitTheTree             ,
  bool pipelineDescendProcessing          ,
  bool pipelineAscendProcessing           ,
  int  smallestProblemSizeForAscendDescend  ,
  int  grainSizeForAscendDescend          ,
  int  smallestProblemSizeForEnterLeaveCell ,
  int  grainSizeForEnterLeaveCell         ,
  int  smallestProblemSizeForTouchFirstLast ,
  int  grainSizeForTouchFirstLast         ,
  int  smallestProblemSizeForSplitLoadStore ,
  int  grainSizeForSplitLoadStore         ,
  int  adapterNumber
):
  _useMulticore(useMultithreading),
  _grainSizeOfUserDefinedRegions(grainSizeOfUserDefinedRegions),
  _adapterNumber(adapterNumber),
  _splitTheTree(splitTheTree),
  _pipelineDescendProcessing(pipelineDescendProcessing),
  _pipelineAscendProcessing(pipelineAscendProcessing),
  _smallestProblemSizeForAscendDescend(smallestProblemSizeForAscendDescend),
  _grainSizeForAscendDescend(grainSizeForAscendDescend),
  _smallestProblemSizeForEnterLeaveCell(smallestProblemSizeForEnterLeaveCell),
  _grainSizeForEnterLeaveCell(grainSizeForEnterLeaveCell),
  _smallestProblemSizeForTouchFirstLast(smallestProblemSizeForTouchFirstLast),
  _grainSizeForTouchFirstLast(grainSizeForTouchFirstLast),
  _smallestProblemSizeForSplitLoadStore(smallestProblemSizeForSplitLoadStore),
  _grainSizeForSplitLoadStore(grainSizeForSplitLoadStore) {
}



peano::datatraversal::autotuning::GrainSize peano::datatraversal::autotuning::OracleForOnePhaseDummy::parallelise(int problemSize, MethodTrace askingMethod) {
  assertion3(problemSize>0,problemSize,peano::datatraversal::autotuning::toString(askingMethod),toString());

  int grainSize           = 0;
  int smallestProblemSize = 0;

  if ( _pipelineAscendProcessing && askingMethod == MethodTrace::PipelineAscendTask ) {
    grainSize           = 1;
    smallestProblemSize = 0;
  }
  else if ( _pipelineDescendProcessing && askingMethod == MethodTrace::PipelineDescendTask ) {
    grainSize           = 1;
    smallestProblemSize = 0;
  }
  else if (
    ( askingMethod==MethodTrace::AscendOnRegularStationaryGrid        ||
      askingMethod==MethodTrace::DescendOnRegularStationaryGrid
    ) &&
    _splitTheTree != SplitVertexReadsOnRegularSubtree::SplitButDoNotParalleliseEvents
  ) {
    grainSize           = _grainSizeForAscendDescend;
    smallestProblemSize = _smallestProblemSizeForAscendDescend;
  }
  else if (
    ( askingMethod==MethodTrace::CallEnterCellOnRegularStationaryGrid ||
      askingMethod==MethodTrace::CallLeaveCellOnRegularStationaryGrid
    ) &&
    _splitTheTree != SplitVertexReadsOnRegularSubtree::SplitButDoNotParalleliseEvents
  ) {
    grainSize           = _grainSizeForEnterLeaveCell;
    smallestProblemSize = _smallestProblemSizeForEnterLeaveCell;
  }
  else if (
    (
      askingMethod==MethodTrace::CallTouchFirstTimeOnRegularStationaryGrid ||
      askingMethod==MethodTrace::CallTouchLastTimeOnRegularStationaryGrid
    ) &&
    _splitTheTree != SplitVertexReadsOnRegularSubtree::SplitButDoNotParalleliseEvents
  ) {
    grainSize           = _grainSizeForTouchFirstLast;
    smallestProblemSize = _smallestProblemSizeForTouchFirstLast;
  }
  else if (
    _splitTheTree != SplitVertexReadsOnRegularSubtree::DoNotSplit  &&
    (
      askingMethod == MethodTrace::SplitLoadVerticesTaskOnRegularStationaryGrid  ||
      askingMethod == MethodTrace::SplitStoreVerticesTaskOnRegularStationaryGrid
    )
  ) {
    grainSize           = _grainSizeForSplitLoadStore;
    smallestProblemSize = _smallestProblemSizeForSplitLoadStore;
  }
  else if (
    askingMethod>=MethodTrace::UserDefined0 && askingMethod<=MethodTrace::UserDefined12
  ) {
    grainSize           = _grainSizeOfUserDefinedRegions;
    smallestProblemSize = _grainSizeOfUserDefinedRegions;
  }


  if (
    _useMulticore
    &&
    problemSize >= smallestProblemSize
  ) {
    assertion4(grainSize<problemSize,grainSize,problemSize,peano::datatraversal::autotuning::toString(askingMethod),toString());
    return GrainSize(grainSize, false, problemSize, askingMethod, this);
  }
  else {
    return GrainSize(0, false, problemSize, askingMethod, this);
  }
}


void peano::datatraversal::autotuning::OracleForOnePhaseDummy::parallelSectionHasTerminated(int problemSize, MethodTrace askingMethod, double costPerProblemElement) {
}


void peano::datatraversal::autotuning::OracleForOnePhaseDummy::loadStatistics(const std::string& filename, int oracleNumber) {
}


void peano::datatraversal::autotuning::OracleForOnePhaseDummy::plotStatistics(std::ostream& out, int oracleNumber) const {
  out << "oracle " << oracleNumber << std::endl
      << toString() << std::endl;
}


peano::datatraversal::autotuning::OracleForOnePhaseDummy::~OracleForOnePhaseDummy() {
}


peano::datatraversal::autotuning::OracleForOnePhase* peano::datatraversal::autotuning::OracleForOnePhaseDummy::createNewOracle() const {
  return new OracleForOnePhaseDummy(
    _useMulticore,
    _grainSizeOfUserDefinedRegions,
    _splitTheTree,
    _pipelineDescendProcessing,
    _pipelineAscendProcessing,
    _smallestProblemSizeForAscendDescend,
    _grainSizeForAscendDescend,
    _smallestProblemSizeForEnterLeaveCell,
    _grainSizeForEnterLeaveCell,
    _smallestProblemSizeForTouchFirstLast,
    _grainSizeForTouchFirstLast,
    _smallestProblemSizeForSplitLoadStore,
    _grainSizeForSplitLoadStore
  );
}


void peano::datatraversal::autotuning::OracleForOnePhaseDummy::deactivateOracle() {
}


void peano::datatraversal::autotuning::OracleForOnePhaseDummy::activateOracle() {
}


std::string peano::datatraversal::autotuning::OracleForOnePhaseDummy::toString(SplitVertexReadsOnRegularSubtree value) {
  switch (value) {
    case SplitVertexReadsOnRegularSubtree::DoNotSplit:
      return "do-not-split";
    case SplitVertexReadsOnRegularSubtree::Split:
      return "split";
    case SplitVertexReadsOnRegularSubtree::SplitButDoNotParalleliseEvents:
      return "split-but-do-not-parallelise-events";
  }
  return "<undef>";
}


std::string peano::datatraversal::autotuning::OracleForOnePhaseDummy::toString() const {
  std::ostringstream msg;

  msg << ",adapter-number="        << _adapterNumber
      << "(multicore="             << _useMulticore
      << ",grain-size-of-user-defined-regions=" << _grainSizeOfUserDefinedRegions
      << ",split-tree="            << toString(_splitTheTree)
      << ",pipeline-descend="      << _pipelineDescendProcessing
      << ",pipeline-ascend="       << _pipelineAscendProcessing
      << ",smallest-problem-size-for-ascend-descend=" << _smallestProblemSizeForAscendDescend
      << ",grain-size-for-ascend-descend=" << _grainSizeForAscendDescend
      << ",smallest-problem-size-for-enter-leave-cell=" << _smallestProblemSizeForEnterLeaveCell
      << ",grain-size-for-enter-leave-cell=" << _grainSizeForEnterLeaveCell
      << ",smallest-problem-size-for-touch-first-last=" << _smallestProblemSizeForTouchFirstLast
      << ",grain-size-for-touch-first-last=" << _grainSizeForTouchFirstLast
      << ",smallest-problem-size-for-split-load-store=" << _smallestProblemSizeForSplitLoadStore
      << ",grain-size-for-split-load-store=" << _grainSizeForSplitLoadStore
      << ")";

  return msg.str();
}
