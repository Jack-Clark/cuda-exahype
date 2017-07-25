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
 
#include "ADERDG2CartesianVTK.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"

#include "exahype/solvers/ADERDGSolver.h"

#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"


#include "kernels/DGBasisFunctions.h"


std::string exahype::plotters::ADERDG2CartesianVerticesVTKAscii::getIdentifier() {
  return "vtk::Cartesian::vertices::ascii";
}


exahype::plotters::ADERDG2CartesianVerticesVTKAscii::ADERDG2CartesianVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
  ADERDG2CartesianVTK(postProcessing,false,false) {
}


std::string exahype::plotters::ADERDG2CartesianVerticesVTKBinary::getIdentifier() {
  return "vtk::Cartesian::vertices::binary";
}


exahype::plotters::ADERDG2CartesianVerticesVTKBinary::ADERDG2CartesianVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2CartesianVTK(postProcessing,true,false) {
}



std::string exahype::plotters::ADERDG2CartesianCellsVTKAscii::getIdentifier() {
  return "vtk::Cartesian::cells::ascii";
}


exahype::plotters::ADERDG2CartesianCellsVTKAscii::ADERDG2CartesianCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2CartesianVTK(postProcessing,false,true) {
}


std::string exahype::plotters::ADERDG2CartesianCellsVTKBinary::getIdentifier() {
 return "vtk::Cartesian::cells::binary";
}


exahype::plotters::ADERDG2CartesianCellsVTKBinary::ADERDG2CartesianCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2CartesianVTK(postProcessing,true,true) {
}


exahype::plotters::ADERDG2CartesianVTK::ADERDG2CartesianVTK(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool plotCells):
  Device(postProcessing),
  _fileCounter(-1),
  _isBinary(isBinary),
  _plotCells(plotCells),
  _order(-1),
  _solverUnknowns(-1),
  _writtenUnknowns(-1),
  _vertexDataWriter(nullptr),
  _cellDataWriter(nullptr),
  _vertexTimeStampDataWriter(nullptr),
  _cellTimeStampDataWriter(nullptr),
  _gridWriter(nullptr),
  _patchWriter(nullptr) {
}


void exahype::plotters::ADERDG2CartesianVTK::init(
  const std::string& filename,
  int                orderPlusOne,
  int                unknowns,
  int                writtenUnknowns,
  const std::string& select
) {
  _filename          = filename;
  _order             = orderPlusOne-1;
  _solverUnknowns    = unknowns;
  _select            = select;
  _patchWriter       = nullptr;
  _writtenUnknowns   = writtenUnknowns;

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


void exahype::plotters::ADERDG2CartesianVTK::startPlotting( double time ) {
  _fileCounter++;

  assertion( _patchWriter==nullptr );

  if (_writtenUnknowns>0) {
    if (_isBinary) {
      _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter());
    }
    else {
      _patchWriter =
        new tarch::plotter::griddata::blockstructured::PatchWriterUnstructured(
          new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter());
    }

    _gridWriter                = _patchWriter->createSinglePatchWriter();
    if (_plotCells) {
      _cellDataWriter          = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);
      _vertexDataWriter        = nullptr;
      _cellTimeStampDataWriter = _patchWriter->createCellDataWriter("time", 1);
    }
    else {
      _cellDataWriter            = nullptr;
      _vertexDataWriter          = _patchWriter->createVertexDataWriter("Q", _writtenUnknowns);
      _vertexTimeStampDataWriter = _patchWriter->createVertexDataWriter("time", 1);
    }

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::ADERDG2CartesianVTK::finishPlotting() {
  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );

    _gridWriter->close();
    if (_vertexDataWriter!=nullptr) _vertexDataWriter->close();
    if (_cellDataWriter!=nullptr)   _cellDataWriter->close();
    if (_vertexTimeStampDataWriter!=nullptr) _vertexTimeStampDataWriter->close();
    if (_cellTimeStampDataWriter!=nullptr)   _cellTimeStampDataWriter->close();

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

  if (_vertexDataWriter!=nullptr)     delete _vertexDataWriter;
  if (_cellDataWriter!=nullptr)       delete _cellDataWriter;
  if (_vertexTimeStampDataWriter!=nullptr)  delete _vertexTimeStampDataWriter;
  if (_cellTimeStampDataWriter!=nullptr)    delete _cellTimeStampDataWriter;
  if (_gridWriter!=nullptr)           delete _gridWriter;
  if (_patchWriter!=nullptr)          delete _patchWriter;

  _vertexDataWriter    = nullptr;
  _cellDataWriter      = nullptr;
  _patchWriter         = nullptr;
  _vertexTimeStampDataWriter = nullptr;
  _cellTimeStampDataWriter   = nullptr;
  _gridWriter          = nullptr;
}



exahype::plotters::ADERDG2CartesianVTK::~ADERDG2CartesianVTK() {
}


void exahype::plotters::ADERDG2CartesianVTK::writeTimeStampDataToPatch( double timeStamp, int vertexIndex, int cellIndex ) {
  if (_writtenUnknowns>0 && _vertexTimeStampDataWriter!=nullptr) {
    dfor(i,_order+1) {
      _vertexTimeStampDataWriter->plotVertex(vertexIndex, timeStamp);
      vertexIndex++;
    }
  }

  if (_writtenUnknowns>0 && _cellTimeStampDataWriter!=nullptr) {
    dfor(i,_order) {
      _cellTimeStampDataWriter->plotCell(cellIndex, timeStamp);
      cellIndex++;
    }
  }
}


void exahype::plotters::ADERDG2CartesianVTK::plotVertexData(
  int firstVertexIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp
) {
  assertion( _vertexDataWriter!=nullptr || _writtenUnknowns==0 );

  double* interpoland = new double[_solverUnknowns];
  double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

  dfor(i,_order+1) {
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = 0.0;
      dfor(ii,_order+1) { // Gauss-Legendre node indices
        int iGauss = peano::utils::dLinearisedWithoutLookup(ii,_order + 1);
        interpoland[unknown] += kernels::equidistantGridProjector1d[_order][ii(1)][i(1)] *
                 kernels::equidistantGridProjector1d[_order][ii(0)][i(0)] *
                 #ifdef Dim3
                 kernels::equidistantGridProjector1d[_order][ii(2)][i(2)] *
                 #endif
                 u[iGauss * _solverUnknowns + unknown];
        assertion3(interpoland[unknown] == interpoland[unknown], offsetOfPatch, sizeOfPatch, iGauss);
      }
    }

    assertion(sizeOfPatch(0)==sizeOfPatch(1));
    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + i.convertScalar<double>()* (sizeOfPatch(0)/(_order)),
      i,
      interpoland,
      value,
      timeStamp
    );

    if (_writtenUnknowns>0) {
      _vertexDataWriter->plotVertex(firstVertexIndex, value, _writtenUnknowns );
    }

    firstVertexIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}


void exahype::plotters::ADERDG2CartesianVTK::plotCellData(
  int firstCellIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp
) {
  assertion( _cellDataWriter!=nullptr || _writtenUnknowns==0 );

  double* interpoland = new double[_solverUnknowns];
  double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

  dfor(i,_order) {
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = kernels::interpolate(
        offsetOfPatch.data(),
        sizeOfPatch.data(),
        (offsetOfPatch + (i.convertScalar<double>()+0.5)* (sizeOfPatch(0)/(_order))).data(),
        _solverUnknowns,
        unknown,
        _order,
        u
      );
    }

    assertion(sizeOfPatch(0)==sizeOfPatch(1));
    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      offsetOfPatch + (i.convertScalar<double>()+0.5)* (sizeOfPatch(0)/(_order)),
      i,
      interpoland,
      value,
      timeStamp
    );

    if (_writtenUnknowns>0) {
      _cellDataWriter->plotCell(firstCellIndex, value, _writtenUnknowns );
    }

    firstCellIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}

void exahype::plotters::ADERDG2CartesianVTK::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2CartesianVTK::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );

    std::pair<int,int> vertexAndCellIndex(0,0);
    if (_writtenUnknowns>0) {
      vertexAndCellIndex = _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order);
    }

    writeTimeStampDataToPatch( timeStamp, vertexAndCellIndex.first, vertexAndCellIndex.second );

    if (_plotCells) {
      plotCellData( vertexAndCellIndex.second, offsetOfPatch, sizeOfPatch, u, timeStamp );
    }
    else {
      plotVertexData( vertexAndCellIndex.first, offsetOfPatch, sizeOfPatch, u, timeStamp );
    }
  }
}
