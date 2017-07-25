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
 
#include "FiniteVolumes2VTKBinary.h"
#include "ADERDG2CartesianVTK.h"

#include "tarch/parallel/Node.h"

#include "exahype/solvers/FiniteVolumesSolver.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"

tarch::logging::Log exahype::plotters::FiniteVolumes2VTKBinary::_log("exahype::plotters::FiniteVolumes2VTKBinary");

exahype::plotters::FiniteVolumes2VTKBinary::FiniteVolumes2VTKBinary(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
  Device(postProcessing),
  _fileCounter(-1),
  _numberOfCellsPerAxis(-1),
  _ghostLayerWidth(ghostLayerWidth),
  _solverUnknowns(-1),
  _writtenUnknowns(-1),
  _patchWriter(nullptr),
  _gridWriter(nullptr),
  _timeStampDataWriter(nullptr),
  _cellDataWriter(nullptr){
}


std::string exahype::plotters::FiniteVolumes2VTKBinary::getIdentifier() {
  return ADERDG2CartesianCellsVTKBinary::getIdentifier();
}


void exahype::plotters::FiniteVolumes2VTKBinary::init(
  const std::string& filename,
  int                numberOfCellsPerAxis,
  int                unknowns,
  int                writtenUnknowns,
  const std::string& select
){
  _filename             = filename;
  _numberOfCellsPerAxis = numberOfCellsPerAxis;
  _solverUnknowns       = unknowns;
  _select               = select;
  _patchWriter          = nullptr;
  _writtenUnknowns      = writtenUnknowns;

  double x;
  x = Parser::getValueFromPropertyString( select, "left" );
  _regionOfInterestLeftBottomFront(0) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
  x = Parser::getValueFromPropertyString( select, "bottom" );
  _regionOfInterestLeftBottomFront(1) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
#ifdef Din3
  x = Parser::getValueFromPropertyString( select, "front" );
  _regionOfInterestLeftBottomFront(2) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
#endif


  x = Parser::getValueFromPropertyString( select, "right" );
  _regionOfInterestRightTopBack(0) = x!=x ? std::numeric_limits<double>::max() : x;
  x = Parser::getValueFromPropertyString( select, "top" );
  _regionOfInterestRightTopBack(1) = x!=x ? std::numeric_limits<double>::max() : x;
#ifdef Dim3
  x = Parser::getValueFromPropertyString( select, "back" );
  _regionOfInterestRightTopBack(2) = x!=x ? std::numeric_limits<double>::max() : x;
#endif
}


void exahype::plotters::FiniteVolumes2VTKBinary::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  if (_writtenUnknowns>0) {
    _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
            new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter());

    _gridWriter          = _patchWriter->createSinglePatchWriter();
    _timeStampDataWriter = _patchWriter->createCellDataWriter("time", 1);
    _cellDataWriter      = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampDataWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::FiniteVolumes2VTKBinary::finishPlotting() {
  _postProcessing->finishPlotting();
  if ( _writtenUnknowns>0 ) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampDataWriter!=nullptr );

    _gridWriter->close();
    _cellDataWriter->close();
    _timeStampDataWriter->close();
    
    std::ostringstream snapshotFileName;
    snapshotFileName << _filename
    #ifdef Parallel
                   << "-rank-" << tarch::parallel::Node::getInstance().getRank()
    #endif
                   << "-" << _fileCounter << ".vtk";

    // See issue #47 for discussion whether to quit program on failure:
    // _patchWriter should raise/throw the C++ Exception or return something in case
    // of failure.
    _patchWriter->writeToFile(snapshotFileName.str());
  }

  if (_cellDataWriter!=nullptr)       delete _cellDataWriter;
  if (_timeStampDataWriter!=nullptr)  delete _timeStampDataWriter;
  if (_gridWriter!=nullptr)           delete _gridWriter;
  if (_patchWriter!=nullptr)          delete _patchWriter;

  _cellDataWriter      = nullptr;
  _patchWriter         = nullptr;
  _timeStampDataWriter = nullptr;
  _gridWriter          = nullptr;
}


exahype::plotters::FiniteVolumes2VTKBinary::~FiniteVolumes2VTKBinary() {
}


void exahype::plotters::FiniteVolumes2VTKBinary::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& cellDescription =
      exahype::solvers::FiniteVolumesSolver::getCellDescription(
          cellDescriptionsIndex,element);

  if (cellDescription.getType()==exahype::solvers::FiniteVolumesSolver::CellDescription::Type::Cell) {
    double* solution = DataHeap::getInstance().getData(cellDescription.getSolution()).data();

    plotPatch(
        cellDescription.getOffset(),
        cellDescription.getSize(), solution,
        cellDescription.getTimeStamp());
  }
}

void exahype::plotters::FiniteVolumes2VTKBinary::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    logDebug("plotPatch(...)","offset of patch: "<<offsetOfPatch
    <<", size of patch: "<<sizeOfPatch
    <<", time stamp: "<<timeStamp);

    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _timeStampDataWriter!=nullptr );

    int cellIndex = _writtenUnknowns==0 ? -1 : _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _numberOfCellsPerAxis).second;

    double* sourceValue = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    dfor(i,_numberOfCellsPerAxis+_ghostLayerWidth) {
      if (tarch::la::allSmaller(i,_numberOfCellsPerAxis+_ghostLayerWidth)
          && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
        if (_writtenUnknowns>0) {
          _timeStampDataWriter->plotCell(cellIndex, timeStamp);
        }

        for (int unknown=0; unknown < _solverUnknowns; unknown++) {
          sourceValue[unknown] =
            u[peano::utils::dLinearisedWithoutLookup(i,_numberOfCellsPerAxis+2*_ghostLayerWidth)*_solverUnknowns+unknown];
        } // !!! Be aware of the "2*_ghostLayerWidth" !!!

        assertion(sizeOfPatch(0)==sizeOfPatch(1));

        _postProcessing->mapQuantities(
          offsetOfPatch,
          sizeOfPatch,
          offsetOfPatch + (i-_ghostLayerWidth).convertScalar<double>()* (sizeOfPatch(0)/(_numberOfCellsPerAxis)),
          i-_ghostLayerWidth,
          sourceValue,
          value,
          timeStamp
        );

        if (_writtenUnknowns>0) {
          _cellDataWriter->plotCell(cellIndex, value, _writtenUnknowns);
        }
        cellIndex++;
      }
    }

    delete[] sourceValue;
    delete[] value;
  }
}
