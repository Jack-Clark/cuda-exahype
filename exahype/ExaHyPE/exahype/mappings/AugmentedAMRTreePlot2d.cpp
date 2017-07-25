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
 
#include "exahype/mappings/AugmentedAMRTreePlot2d.h"

#include <sstream>

#include "peano/grid/CellFlags.h"
#include "peano/utils/Loop.h"

#include "tarch/la/VectorScalarOperations.h"

#include "exahype/solvers/Solver.h"

#ifdef Parallel
#include "tarch/parallel/Node.h"
#endif

tarch::logging::Log exahype::mappings::AugmentedAMRTreePlot2d::_log(
    "exahype::mappings::AugmentedAMRTreePlot2d");

int exahype::mappings::AugmentedAMRTreePlot2d::_snapshotCounter        = 0;
double exahype::mappings::AugmentedAMRTreePlot2d::SqueezeZAxis         = 4.0;
double exahype::mappings::AugmentedAMRTreePlot2d::TreeConnectionsValue = -100.0;

peano::CommunicationSpecification
exahype::mappings::AugmentedAMRTreePlot2d::communicationSpecification() {
  return peano::CommunicationSpecification::getPessimisticSpecification(true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::touchVertexLastTimeSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::touchVertexFirstTimeSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial,true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::enterCellSpecification() {
  return peano::MappingSpecification(peano::MappingSpecification::WholeTree,
                                     peano::MappingSpecification::Serial,true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::leaveCellSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::AvoidFineGridRaces,true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::ascendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

peano::MappingSpecification
exahype::mappings::AugmentedAMRTreePlot2d::descendSpecification() {
  return peano::MappingSpecification(
      peano::MappingSpecification::Nop,
      peano::MappingSpecification::RunConcurrentlyOnFineGrid,true);
}

std::map<tarch::la::Vector<DIMENSIONS + 1, double>, int,
         tarch::la::VectorCompare<DIMENSIONS + 1> >
    exahype::mappings::AugmentedAMRTreePlot2d::_vertex2IndexMap;
std::map<tarch::la::Vector<DIMENSIONS + 1, double>, int,
         tarch::la::VectorCompare<DIMENSIONS + 1> >
    exahype::mappings::AugmentedAMRTreePlot2d::_cellCenter2IndexMap;

std::map<int,double> exahype::mappings::AugmentedAMRTreePlot2d::_level2OffsetMap;

exahype::mappings::AugmentedAMRTreePlot2d::AugmentedAMRTreePlot2d()
    : _vtkWriter(0),
      _vertexWriter(0),
      _cellWriter(0),
      _cellNumberWriter(0),
      _cellTypeWriter(0),
      _cellDescriptionIndexWriter(0),
      _cellRefinementEventWriter(0),
      _cellDataWriter(0),
      _limiterStatusWriter(0),
      _cellCounter(0) {}

exahype::mappings::AugmentedAMRTreePlot2d::~AugmentedAMRTreePlot2d() {}

#if defined(SharedMemoryParallelisation)
exahype::mappings::AugmentedAMRTreePlot2d::AugmentedAMRTreePlot2d(
    const AugmentedAMRTreePlot2d& masterThread)
    : _vtkWriter(masterThread._vtkWriter),
      _vertexWriter(masterThread._vertexWriter),
      _cellWriter(masterThread._cellWriter),
      _cellNumberWriter(masterThread._cellNumberWriter),
      _cellTypeWriter(masterThread._cellTypeWriter),
      _cellDescriptionIndexWriter(masterThread._cellDescriptionIndexWriter),
      _cellRefinementEventWriter(masterThread._cellRefinementEventWriter),
      _cellDataWriter(masterThread._cellDataWriter),
      _limiterStatusWriter(masterThread._limiterStatusWriter),
      _cellCounter(0) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithWorkerThread(
    const AugmentedAMRTreePlot2d& workerThread) {}
#endif

void exahype::mappings::AugmentedAMRTreePlot2d::plotVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x, int level) {
  tarch::la::Vector<DIMENSIONS + 1, double> y;
  for (int i = 0; i < DIMENSIONS; i++) y(i) = x(i);
  y(DIMENSIONS) = level;

  tarch::la::Vector<DIMENSIONS + 1, double> plotY = y;
  plotY(DIMENSIONS) = level / SqueezeZAxis;

#if DIMENSIONS == 2
  if (_vertex2IndexMap.find(y) == _vertex2IndexMap.end()) {
    _vertex2IndexMap[y] = _vertexWriter->plotVertex(plotY);
  }
#else
  logError("plotVertex",
           "This mapping can only be used for two-dimensional problems "
           "(DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::createHangingVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
#if DIMENSIONS == 2
  double offsetZ = coarseGridVerticesEnumerator.getLevel()+1;
  plotVertex(fineGridVertex, fineGridX, offsetZ);
#else
  logError("createHangingVertex",
           "This mapping can only be used for two-dimensional problems "
           "(DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::destroyHangingVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::createInnerVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::createBoundaryVertex(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::destroyVertex(
    const exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::createCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::destroyCell(
    const exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

#ifdef Parallel
void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithNeighbour(
    exahype::Vertex& vertex, const exahype::Vertex& neighbour, int fromRank,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::prepareSendToNeighbour(
    exahype::Vertex& vertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::prepareCopyToRemoteNode(
    exahype::Vertex& localVertex, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::prepareCopyToRemoteNode(
    exahype::Cell& localCell, int toRank,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Vertex& localVertex,
        const exahype::Vertex& masterOrWorkerVertex, int fromRank,
        const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {
}

void exahype::mappings::AugmentedAMRTreePlot2d::
    mergeWithRemoteDataDueToForkOrJoin(
        exahype::Cell& localCell, const exahype::Cell& masterOrWorkerCell,
        int fromRank, const tarch::la::Vector<DIMENSIONS, double>& x,
        const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}

bool exahype::mappings::AugmentedAMRTreePlot2d::prepareSendToWorker(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell,
    int worker) {
  return false;
}

void exahype::mappings::AugmentedAMRTreePlot2d::prepareSendToMaster(
    exahype::Cell& localCell, exahype::Vertex* vertices,
    const peano::grid::VertexEnumerator& verticesEnumerator,
    const exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    const exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithMaster(
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
    exahype::State& masterState) {}

void exahype::mappings::AugmentedAMRTreePlot2d::receiveDataFromMaster(
    exahype::Cell& receivedCell, exahype::Vertex* receivedVertices,
    const peano::grid::VertexEnumerator& receivedVerticesEnumerator,
    exahype::Vertex* const receivedCoarseGridVertices,
    const peano::grid::VertexEnumerator& receivedCoarseGridVerticesEnumerator,
    exahype::Cell& receivedCoarseGridCell,
    exahype::Vertex* const workersCoarseGridVertices,
    const peano::grid::VertexEnumerator& workersCoarseGridVerticesEnumerator,
    exahype::Cell& workersCoarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithWorker(
    exahype::Cell& localCell, const exahype::Cell& receivedMasterCell,
    const tarch::la::Vector<DIMENSIONS, double>& cellCentre,
    const tarch::la::Vector<DIMENSIONS, double>& cellSize, int level) {}

void exahype::mappings::AugmentedAMRTreePlot2d::mergeWithWorker(
    exahype::Vertex& localVertex, const exahype::Vertex& receivedMasterVertex,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, double>& h, int level) {}
#endif

void exahype::mappings::AugmentedAMRTreePlot2d::touchVertexFirstTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {
#if DIMENSIONS == 2
  double offsetZ = coarseGridVerticesEnumerator.getLevel()+1;
  plotVertex(fineGridVertex, fineGridX, offsetZ);
#else
  logError("touchVertexFirstTime",
           "This mapping can only be used for two-dimensional problems "
           "(DIMENSIONS==2).")
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::touchVertexLastTime(
    exahype::Vertex& fineGridVertex,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridX,
    const tarch::la::Vector<DIMENSIONS, double>& fineGridH,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfVertex) {}

void exahype::mappings::AugmentedAMRTreePlot2d::enterCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {
#if DIMENSIONS == 2
    double offsetZ = coarseGridVerticesEnumerator.getLevel()+1;

    int vertexIndex[TWO_POWER_D];
    tarch::la::Vector<DIMENSIONS + 1, double> currentVertexPosition;
    currentVertexPosition(DIMENSIONS) = offsetZ;

    dfor2(i) for (int d = 0; d < DIMENSIONS; d++) {
      currentVertexPosition(d) = fineGridVerticesEnumerator.getVertexPosition(i)(d);
    }
    assertion2(
        _vertex2IndexMap.find(currentVertexPosition) != _vertex2IndexMap.end(),
        currentVertexPosition,
        fineGridVertices[fineGridVerticesEnumerator(i)].toString());
    vertexIndex[iScalar] = _vertex2IndexMap[currentVertexPosition];
    enddforx

    _cellWriter->plotQuadrangle(vertexIndex);

    const int cellIndex = _cellWriter->plotQuadrangle(vertexIndex);

    _cellNumberWriter->plotCell(cellIndex, _cellCounter);

    _cellDescriptionIndexWriter->plotCell(
        cellIndex,
        static_cast<int>(fineGridCell.getCellDescriptionsIndex()));

    if (fineGridCell.isInitialised() &&
        exahype::solvers::ADERDGSolver::Heap::getInstance().getData(fineGridCell.getCellDescriptionsIndex()).size() > 0) {
      int solverNumber = 0;
      bool solverFound = false;

      for (auto& pFine : exahype::solvers::ADERDGSolver::Heap::getInstance().getData(
               fineGridCell.getCellDescriptionsIndex())) {
        if (pFine.getSolverNumber() == solverNumber) {
          _cellTypeWriter->plotCell(cellIndex,
                                    static_cast<int>(pFine.getType()));
          _cellRefinementEventWriter->plotCell(
              cellIndex, static_cast<int>(pFine.getRefinementEvent()));
          _cellDataWriter->plotCell(
              cellIndex,
              2 * static_cast<int>(pFine.getSolution() > -1) +
                  static_cast<int>(pFine.getExtrapolatedPredictor() > -1));
          _limiterStatusWriter->plotCell(
              cellIndex, static_cast<int>(pFine.getLimiterStatus()));
          solverFound = true;
        }
      }

      if (!solverFound) {
        _cellTypeWriter->plotCell(cellIndex, -1);
        _cellRefinementEventWriter->plotCell(cellIndex, -1);
        _cellDataWriter->plotCell(cellIndex, 0);
        _limiterStatusWriter->plotCell(cellIndex, -1);
      }

    } else {
      _cellTypeWriter->plotCell(
          cellIndex,
          static_cast<int>(fineGridCell.getCellDescriptionsIndex()));
      _cellRefinementEventWriter->plotCell(cellIndex, -1);
      _cellDataWriter->plotCell(cellIndex, 0);
      _limiterStatusWriter->plotCell(cellIndex, -1);
    }

    _cellCounter++;
#endif
}

void exahype::mappings::AugmentedAMRTreePlot2d::leaveCell(
    exahype::Cell& fineGridCell, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell,
    const tarch::la::Vector<DIMENSIONS, int>& fineGridPositionOfCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::beginIteration(
    exahype::State& solverState) {
  assertion(_vtkWriter == 0);

  _cellCounter = 0;

  _vtkWriter = new UsedWriter();

  _vertexWriter = _vtkWriter->createVertexWriter();
  _cellWriter   = _vtkWriter->createCellWriter();

  _cellNumberWriter = _vtkWriter->createCellDataWriter("cell-number", 1);
  _cellTypeWriter   = _vtkWriter->createCellDataWriter(
      "cell-type(NoPatch=-1,Erased=0,Ancestor=1,EmptyAncestor=2,Cell=3,"
      "Descendant=4,EmptyDescendant=5)",
      1);
  _cellDescriptionIndexWriter =
      _vtkWriter->createCellDataWriter("NoPatch=-1,ValidPatch>=0", 1);
  _cellRefinementEventWriter = _vtkWriter->createCellDataWriter(
      "refinement-event(NoPatch=-1,None=0,ErasingChildrenReq=1,"
      "ErasingChildren=2,ChangeToDescendantsReq=3,"
      "ChangeToDescendants=4,RefReq=5,Ref=6,"
      "DeaugChildrenReq=7,DeaugChildren=8,"
      "AugReq=9,Aug=10)"
      ,1);
  _cellDataWriter = _vtkWriter->createCellDataWriter(
      "Data-on-Patch(None=0,OnlyFaceData=1,VolumeAndFaceData=3)", 1);
  _limiterStatusWriter = _vtkWriter->createCellDataWriter(
      "Limiter-Status(Ok=0,NeighbourOfNeighbourOfTroubled=1,NeighbourOfTroubled=2,Troubled=3)", 1);
}

void exahype::mappings::AugmentedAMRTreePlot2d::endIteration(
    exahype::State& solverState) {
  _vertexWriter->close();
  _cellWriter->close();

  _cellTypeWriter->close();
  _cellDescriptionIndexWriter->close();
  _cellRefinementEventWriter->close();
  _cellDataWriter->close();
  _limiterStatusWriter->close();
  _cellNumberWriter->close();

  delete _vertexWriter;
  delete _cellWriter;

  delete _cellTypeWriter;
  delete _cellDescriptionIndexWriter;
  delete _cellRefinementEventWriter;
  delete _cellNumberWriter;
  delete _cellDataWriter;
  delete _limiterStatusWriter;

  _vertexWriter = 0;
  _cellWriter = 0;

  _cellNumberWriter = 0;
  _cellTypeWriter = 0;
  _cellDescriptionIndexWriter = 0;
  _cellRefinementEventWriter = 0;
  _cellDataWriter = 0;

  std::ostringstream snapshotFileName;
  snapshotFileName << "tree"
#ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
#endif
                   << "-" << _snapshotCounter << ".vtk";
  _vtkWriter->writeToFile(snapshotFileName.str());

  _snapshotCounter++;

  _vertex2IndexMap.clear();

  delete _vtkWriter;
  _vtkWriter = 0;
}

void exahype::mappings::AugmentedAMRTreePlot2d::descend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}

void exahype::mappings::AugmentedAMRTreePlot2d::ascend(
    exahype::Cell* const fineGridCells, exahype::Vertex* const fineGridVertices,
    const peano::grid::VertexEnumerator& fineGridVerticesEnumerator,
    exahype::Vertex* const coarseGridVertices,
    const peano::grid::VertexEnumerator& coarseGridVerticesEnumerator,
    exahype::Cell& coarseGridCell) {}
