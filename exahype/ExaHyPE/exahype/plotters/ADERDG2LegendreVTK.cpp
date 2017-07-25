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
 
#include "ADERDG2LegendreVTK.h"
#include "tarch/parallel/Node.h"

#include "kernels/DGMatrices.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGBasisFunctions.h"

#include "peano/utils/Loop.h"

#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"

#include "exahype/solvers/ADERDGSolver.h"




std::string exahype::plotters::ADERDG2LegendreVerticesVTKAscii::getIdentifier() {
  return "vtk::Legendre::vertices::ascii";
}


exahype::plotters::ADERDG2LegendreVerticesVTKAscii::ADERDG2LegendreVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2LegendreVTK(postProcessing,false,false) {
}


std::string exahype::plotters::ADERDG2LegendreVerticesVTKBinary::getIdentifier() {
  return "vtk::Legendre::vertices::binary";
}


exahype::plotters::ADERDG2LegendreVerticesVTKBinary::ADERDG2LegendreVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2LegendreVTK(postProcessing,true,false) {
}



std::string exahype::plotters::ADERDG2LegendreCellsVTKAscii::getIdentifier() {
  return "vtk::Legendre::cells::ascii";
}


exahype::plotters::ADERDG2LegendreCellsVTKAscii::ADERDG2LegendreCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2LegendreVTK(postProcessing,false,true) {
}


std::string exahype::plotters::ADERDG2LegendreCellsVTKBinary::getIdentifier() {
 return "vtk::Legendre::cells::binary";
}


exahype::plotters::ADERDG2LegendreCellsVTKBinary::ADERDG2LegendreCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2LegendreVTK(postProcessing,true,true) {
}



exahype::plotters::ADERDG2LegendreVTK::ADERDG2LegendreVTK(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool plotCells):
  Device(postProcessing),
  _fileCounter(-1),
  _isBinary(isBinary),
  _plotCells(plotCells),
  _gridWriter(nullptr),
  _vertexWriter(nullptr),
  _cellWriter(nullptr),
  _vertexTimeStampDataWriter(nullptr),
  _cellTimeStampDataWriter(nullptr),
  _vertexDataWriter(nullptr),
  _cellDataWriter(nullptr) {
}


void exahype::plotters::ADERDG2LegendreVTK::init(
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
  _writtenUnknowns   = writtenUnknowns;

  double x;
  x = Parser::getValueFromPropertyString( select, "left" );
  _regionOfInterestLeftBottomFront(0) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
  x = Parser::getValueFromPropertyString( select, "bottom" );
  _regionOfInterestLeftBottomFront(1) = x!=x ? -std::numeric_limits<double>::max() : x; // "-", min
#ifdef Dim3
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


void exahype::plotters::ADERDG2LegendreVTK::startPlotting( double time ) {
  _fileCounter++;

  if (_writtenUnknowns>0) {
    if (_isBinary) {
      _gridWriter = new tarch::plotter::griddata::unstructured::vtk::VTKBinaryFileWriter();
    }
    else {
      _gridWriter = new tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter();
    }

    _vertexWriter                = _gridWriter->createVertexWriter();
    _cellWriter                  = _gridWriter->createCellWriter();

    if (_plotCells) {
      _cellDataWriter          = _gridWriter->createCellDataWriter("Q", _writtenUnknowns);
      _vertexDataWriter        = nullptr;
      _cellTimeStampDataWriter = _gridWriter->createCellDataWriter("time", 1);
    }
    else {
      _cellDataWriter          = nullptr;
      _vertexDataWriter        = _gridWriter->createVertexDataWriter("Q", _writtenUnknowns);
      _vertexTimeStampDataWriter = _gridWriter->createVertexDataWriter("time", 1);
    }


    assertion( _gridWriter!=nullptr );
    assertion( _vertexWriter!=nullptr );
    assertion( _cellWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::ADERDG2LegendreVTK::finishPlotting() {
  _postProcessing->finishPlotting();

  if ( _writtenUnknowns>0 ) {
    assertion( _gridWriter!=nullptr );

    _vertexWriter->close();
    _cellWriter->close();
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
    _gridWriter->writeToFile(snapshotFileName.str());
  }

  if (_vertexDataWriter!=nullptr)    delete _vertexDataWriter;
  if (_cellDataWriter!=nullptr)      delete _cellDataWriter;
  if (_vertexWriter!=nullptr)        delete _vertexWriter;
  if (_cellWriter!=nullptr)          delete _cellWriter;
  if (_vertexTimeStampDataWriter!=nullptr) delete _vertexTimeStampDataWriter;
  if (_cellTimeStampDataWriter!=nullptr)   delete _cellTimeStampDataWriter;
  if (_gridWriter!=nullptr)          delete _gridWriter;

  _vertexDataWriter     = nullptr;
  _cellDataWriter       = nullptr;
  _vertexWriter         = nullptr;
  _cellWriter           = nullptr;
  _vertexTimeStampDataWriter  = nullptr;
  _cellTimeStampDataWriter    = nullptr;
  _gridWriter           = nullptr;
}


exahype::plotters::ADERDG2LegendreVTK::~ADERDG2LegendreVTK() {
}


void exahype::plotters::ADERDG2LegendreVTK::writeTimeStampDataToPatch( double timeStamp, int vertexIndex, int cellIndex ) {
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


std::pair<int,int> exahype::plotters::ADERDG2LegendreVTK::plotLegendrePatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch
) {
  int firstVertex = -1;
  int firstCell   = -1;

  if (_writtenUnknowns>0) {
    assertion(_vertexWriter!=nullptr);
    dfor(i,_order+1) {
      tarch::la::Vector<2, double> p;

      //p = offsetOfPatch + tarch::la::multiplyComponents( i.convertScalar<double>(), sizeOfPatch) * (1.0/_order);

      for (int d=0; d<DIMENSIONS; d++) {
        p(d) = offsetOfPatch(d) + kernels::gaussLegendreNodes[_order][i(d)] * sizeOfPatch(d);
      }

      const int newVertexNumber = _vertexWriter->plotVertex(p);
      firstVertex = firstVertex==-1 ? newVertexNumber : firstVertex;
    }

    assertion(_cellWriter!=nullptr);
    dfor(i,_order) {
      #ifdef Dim2
      int cellsVertexIndices[4];
      cellsVertexIndices[0] = firstVertex + (i(0)+0) + (i(1)+0) * (_order+1);
      cellsVertexIndices[1] = firstVertex + (i(0)+1) + (i(1)+0) * (_order+1);
      cellsVertexIndices[2] = firstVertex + (i(0)+0) + (i(1)+1) * (_order+1);
      cellsVertexIndices[3] = firstVertex + (i(0)+1) + (i(1)+1) * (_order+1);
      const int newCellNumber = _cellWriter->plotQuadrangle(cellsVertexIndices);
      firstCell = firstCell==-1 ? newCellNumber : firstCell;
      #elif Dim3
      assertionMsg( false, "not implemented yet" );
      int cellsVertexIndices[8];
      cellsVertexIndices[0] = firstVertex + (i(0)+0) + (i(1)+0) * (_order+1) + (i(2)+0) * (_order+1) * (_order+1);
      cellsVertexIndices[1] = firstVertex + (i(0)+1) + (i(1)+0) * (_order+1) + (i(2)+0) * (_order+1) * (_order+1);
      cellsVertexIndices[2] = firstVertex + (i(0)+0) + (i(1)+1) * (_order+1) + (i(2)+0) * (_order+1) * (_order+1);
      cellsVertexIndices[3] = firstVertex + (i(0)+1) + (i(1)+1) * (_order+1) + (i(2)+0) * (_order+1) * (_order+1);
      cellsVertexIndices[4] = firstVertex + (i(0)+0) + (i(1)+0) * (_order+1) + (i(2)+1) * (_order+1) * (_order+1);
      cellsVertexIndices[5] = firstVertex + (i(0)+1) + (i(1)+0) * (_order+1) + (i(2)+1) * (_order+1) * (_order+1);
      cellsVertexIndices[6] = firstVertex + (i(0)+0) + (i(1)+1) * (_order+1) + (i(2)+1) * (_order+1) * (_order+1);
      cellsVertexIndices[7] = firstVertex + (i(0)+1) + (i(1)+1) * (_order+1) + (i(2)+1) * (_order+1) * (_order+1);
      const int newCellNumber = _cellWriter->plotHexahedron(cellsVertexIndices);
      firstCell = firstCell==-1 ? newCellNumber : firstCell;
      #endif
    }
  }

  return std::pair<int,int>(firstVertex,firstCell);
}


void exahype::plotters::ADERDG2LegendreVTK::plotVertexData(
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
    // This is inefficient but works. We could look it up directly from the arrays
    tarch::la::Vector<DIMENSIONS, double> p;
    for (int d=0; d<DIMENSIONS; d++) {
      p(d) = offsetOfPatch(d) + kernels::gaussLegendreNodes[_order][i(d)] * sizeOfPatch(d);
    }
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = kernels::interpolate(
        offsetOfPatch.data(),
        sizeOfPatch.data(),
        p.data(), // das ist die Position
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
      p,
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


void exahype::plotters::ADERDG2LegendreVTK::plotCellData(
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
    // This is inefficient but works. We could look it up directly from the arrays
    tarch::la::Vector<DIMENSIONS, double> p;
    for (int d=0; d<DIMENSIONS; d++) {
      p(d) = offsetOfPatch(d) + (kernels::gaussLegendreNodes[_order][i(d)]+kernels::gaussLegendreNodes[_order][i(d)+1]) * sizeOfPatch(d)/2.0;
    }
    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = kernels::interpolate(
        offsetOfPatch.data(),
        sizeOfPatch.data(),
        p.data(), // das ist die Position
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
      p,
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

void exahype::plotters::ADERDG2LegendreVTK::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2LegendreVTK::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    assertion( _writtenUnknowns==0 || _vertexWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _cellWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _vertexTimeStampDataWriter!=nullptr );

    std::pair<int,int> vertexAndCellIndex = plotLegendrePatch(offsetOfPatch, sizeOfPatch);

    writeTimeStampDataToPatch( timeStamp, vertexAndCellIndex.first, vertexAndCellIndex.second );

    if (_plotCells) {
      plotCellData( vertexAndCellIndex.second, offsetOfPatch, sizeOfPatch, u, timeStamp );
    }
    else {
      plotVertexData( vertexAndCellIndex.first, offsetOfPatch, sizeOfPatch, u, timeStamp );
    }
  }
}
