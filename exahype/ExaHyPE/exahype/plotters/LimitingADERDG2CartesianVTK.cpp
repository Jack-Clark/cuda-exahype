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
 
#include "LimitingADERDG2CartesianVTK.h"
#include "tarch/parallel/Node.h"

// @todo 16/05/03:Dominic Etienne Charreir Plotter depends now on kernels.
// Should thus be placed in kernel module or the solver
// should provide a function that computes solution values
// at equidistant grid points
#include "kernels/DGMatrices.h"
#include "peano/utils/Loop.h"


#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKBinaryFileWriter.h"


#include "kernels/DGBasisFunctions.h"

#include "exahype/solvers/LimitingADERDGSolver.h"


std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTKAscii::getIdentifier() {
  return "vtk::Cartesian::vertices::limited::ascii";
}


exahype::plotters::LimitingADERDG2CartesianVerticesVTKAscii::LimitingADERDG2CartesianVerticesVTKAscii(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
  LimitingADERDG2CartesianVTK(postProcessing,ghostLayerWidth,false,false) {
}


std::string exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinary::getIdentifier() {
  return "vtk::Cartesian::vertices::limited::binary";
}


exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinary::LimitingADERDG2CartesianVerticesVTKBinary(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTK(postProcessing,ghostLayerWidth,true,false) {
}



std::string exahype::plotters::LimitingADERDG2CartesianCellsVTKAscii::getIdentifier() {
  return "vtk::Cartesian::cells::limited::ascii";
}


exahype::plotters::LimitingADERDG2CartesianCellsVTKAscii::LimitingADERDG2CartesianCellsVTKAscii(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTK(postProcessing,ghostLayerWidth,false,true) {
}


std::string exahype::plotters::LimitingADERDG2CartesianCellsVTKBinary::getIdentifier() {
 return "vtk::Cartesian::cells::limited::binary";
}


exahype::plotters::LimitingADERDG2CartesianCellsVTKBinary::LimitingADERDG2CartesianCellsVTKBinary(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth):
    LimitingADERDG2CartesianVTK(postProcessing,ghostLayerWidth,true,true) {
}

tarch::logging::Log exahype::plotters::LimitingADERDG2CartesianVTK::_log("exahype::plotters::LimitingADERDG2CartesianVTK");

exahype::plotters::LimitingADERDG2CartesianVTK::LimitingADERDG2CartesianVTK(
    exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
    const int ghostLayerWidth,
    const bool isBinary, const bool plotCells)
  :
  Device(postProcessing),
  _fileCounter(-1),
  _isBinary(isBinary),
  _plotCells(plotCells),
  _order(-1),
  _solverUnknowns(-1),
  _writtenUnknowns(-1),
  _ghostLayerWidth(ghostLayerWidth),
  _gridWriter(nullptr),
  _patchWriter(nullptr),
  _vertexDataWriter(nullptr),
  _cellDataWriter(nullptr),
  _timeStampVertexDataWriter(nullptr),
  _timeStampCellDataWriter(nullptr),
  _cellLimiterStatusWriter(nullptr),
  _vertexLimiterStatusWriter(nullptr){
}


void exahype::plotters::LimitingADERDG2CartesianVTK::init(
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


void exahype::plotters::LimitingADERDG2CartesianVTK::startPlotting( double time ) {
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

    _gridWriter                  = _patchWriter->createSinglePatchWriter();
    if (_plotCells) {
      _cellDataWriter            = _patchWriter->createCellDataWriter("Q", _writtenUnknowns);
      _vertexDataWriter          = nullptr;

      _cellLimiterStatusWriter   = _patchWriter->createCellDataWriter("Limiter-Status(0-O,1-NNT,2-NT,3-T)", 1);
      _vertexLimiterStatusWriter = nullptr;
    }
    else {
      _cellDataWriter            = nullptr;
      _vertexDataWriter          = _patchWriter->createVertexDataWriter("Q", _writtenUnknowns);

      _cellLimiterStatusWriter   = nullptr;
      _vertexLimiterStatusWriter = _patchWriter->createVertexDataWriter("Limiter-Status(0-O,1-NNT,2-NT,3-T)", 1);
    }
    _timeStampVertexDataWriter = _patchWriter->createVertexDataWriter("time", 1);
//    _timeStampCellDataWriter   = _patchWriter->createCellDataWriter("time", 1);

    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampVertexDataWriter!=nullptr );
//    assertion( _timeStampCellDataWriter!=nullptr );
  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::LimitingADERDG2CartesianVTK::finishPlotting() {
  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {
    assertion( _patchWriter!=nullptr );
    assertion( _gridWriter!=nullptr );
    assertion( _timeStampVertexDataWriter!=nullptr );

    _gridWriter->close();
//    if (_timeStampCellDataWriter!=nullptr) _timeStampCellDataWriter->close();
    if (_vertexDataWriter!=nullptr)        _vertexDataWriter->close();
    if (_cellDataWriter!=nullptr)          _cellDataWriter->close();
    if (_cellLimiterStatusWriter!=nullptr) _cellLimiterStatusWriter->close();
    if (_vertexLimiterStatusWriter!=nullptr) _vertexLimiterStatusWriter->close();
    _timeStampVertexDataWriter->close();

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

  if (_vertexDataWriter!=nullptr)          delete _vertexDataWriter;
  if (_cellDataWriter!=nullptr)            delete _cellDataWriter;
  if (_timeStampVertexDataWriter!=nullptr) delete _timeStampVertexDataWriter;
  if (_timeStampCellDataWriter!=nullptr)   delete _timeStampCellDataWriter;
  if (_cellLimiterStatusWriter!=nullptr)   delete _cellLimiterStatusWriter;
  if (_vertexLimiterStatusWriter!=nullptr) delete _vertexLimiterStatusWriter;
  if (_gridWriter!=nullptr)                delete _gridWriter;
  if (_patchWriter!=nullptr)               delete _patchWriter;

  _vertexDataWriter          = nullptr;
  _cellDataWriter            = nullptr;
  _patchWriter               = nullptr;
  _timeStampVertexDataWriter = nullptr;
  _timeStampCellDataWriter   = nullptr;
  _cellLimiterStatusWriter   = nullptr;
  _vertexLimiterStatusWriter = nullptr;
  _gridWriter                = nullptr;
}



exahype::plotters::LimitingADERDG2CartesianVTK::~LimitingADERDG2CartesianVTK() {
}


void exahype::plotters::LimitingADERDG2CartesianVTK::writeTimeStampDataToADERDGPatch( double timeStamp, int vertexIndex ) {
  if (_writtenUnknowns>0) {
    dfor(i,_order+1) {
      _timeStampVertexDataWriter->plotVertex(vertexIndex, timeStamp);
      vertexIndex++;
    }
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTK::plotVertexData(
  int firstVertexIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp,
  const int limiterStatusAsInt
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

    _vertexLimiterStatusWriter->plotVertex(firstVertexIndex, static_cast<double>(limiterStatusAsInt));

    firstVertexIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}


void exahype::plotters::LimitingADERDG2CartesianVTK::plotCellData(
  int firstCellIndex,
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
  double* u,
  double timeStamp,
  const int limiterStatusAsInt
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

    _cellLimiterStatusWriter->plotCell(firstCellIndex, static_cast<double>(limiterStatusAsInt));

    firstCellIndex++;
  }

  if (interpoland!=nullptr)  delete[] interpoland;
  if (value!=nullptr)        delete[] value;
}

void exahype::plotters::LimitingADERDG2CartesianVTK::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& solverPatch = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (solverPatch.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    switch(solverPatch.getLimiterStatus()) {
      case exahype::records::ADERDGCellDescription::LimiterStatus::Troubled:
      case exahype::records::ADERDGCellDescription::LimiterStatus::NeighbourIsTroubledCell: // TODO(Dominic): Plot FVM solution instead
      case exahype::records::ADERDGCellDescription::LimiterStatus::Ok:
      case exahype::records::ADERDGCellDescription::LimiterStatus::NeighbourIsNeighbourOfTroubledCell: {
        double* solverSolution = DataHeap::getInstance().getData(solverPatch.getSolution()).data();

        plotADERDGPatch(
            solverPatch.getOffset(),
            solverPatch.getSize(), solverSolution,
            solverPatch.getCorrectorTimeStamp(),
            static_cast<int>(solverPatch.getLimiterStatus()));
      } break;
//      case exahype::records::ADERDGCellDescription::LimiterStatus::Troubled:
//      case exahype::records::ADERDGCellDescription::LimiterStatus::NeighbourIsTroubledCell: {
//        auto* limitingADERDGSolver =
//            static_cast<exahype::solvers::LimitingADERDGSolver*>(
//                exahype::solvers::RegisteredSolvers[solverPatch.getSolverNumber()]);
//
//        const int limiterElement =
//            limitingADERDGSolver->tryGetLimiterElement(cellDescriptionsIndex,solverPatch.getSolverNumber());
//        auto& limiterPatch =
//            exahype::solvers::FiniteVolumesSolver::getCellDescription(cellDescriptionsIndex,limiterElement);
//
//        double* limiterSolution = DataHeap::getInstance().getData(limiterPatch.getSolution()).data();
//
//        plotFiniteVolumesPatch(
//            limiterPatch.getOffset(),
//            limiterPatch.getSize(), limiterSolution,
//            limiterPatch.getTimeStamp());
//      } break;
    }
  }
}


void exahype::plotters::LimitingADERDG2CartesianVTK::plotADERDGPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp,
    const int limiterStatusAsInt) {
  if (
    tarch::la::allSmaller(_regionOfInterestLeftBottomFront,offsetOfPatch+sizeOfPatch)
    &&
    tarch::la::allGreater(_regionOfInterestRightTopBack,offsetOfPatch)
  ) {
    assertion( _writtenUnknowns==0 || _patchWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _gridWriter!=nullptr );
    assertion( _writtenUnknowns==0 || _timeStampVertexDataWriter!=nullptr );

    std::pair<int,int> vertexAndCellIndex(0,0);
    if (_writtenUnknowns>0) {
      vertexAndCellIndex = _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, _order);
    }

    writeTimeStampDataToADERDGPatch( timeStamp, vertexAndCellIndex.first );

    if (_plotCells) {
      plotCellData( vertexAndCellIndex.second, offsetOfPatch, sizeOfPatch, u, timeStamp, limiterStatusAsInt );
    }
    else {
      plotVertexData( vertexAndCellIndex.first, offsetOfPatch, sizeOfPatch, u, timeStamp, limiterStatusAsInt );
    }
  }
}

void exahype::plotters::LimitingADERDG2CartesianVTK::plotFiniteVolumesPatch(
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
    assertion( _writtenUnknowns==0 || _timeStampVertexDataWriter!=nullptr );

    const int numberOfCellsPerAxis = 2*_order+1;

    int cellIndex = _writtenUnknowns==0 ? -1 : _gridWriter->plotPatch(offsetOfPatch, sizeOfPatch, numberOfCellsPerAxis).second;

    double* sourceValue = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    dfor(i,numberOfCellsPerAxis+_ghostLayerWidth) {
      if (tarch::la::allSmaller(i,numberOfCellsPerAxis+_ghostLayerWidth)
          && tarch::la::allGreater(i,_ghostLayerWidth-1)) {
        if (_writtenUnknowns>0) {
          _timeStampCellDataWriter->plotCell(cellIndex, timeStamp);
        }

        for (int unknown=0; unknown < _solverUnknowns; unknown++) {
          sourceValue[unknown] =
            u[peano::utils::dLinearisedWithoutLookup(i,numberOfCellsPerAxis+2*_ghostLayerWidth)*_solverUnknowns+unknown];
        } // !!! Be aware of the "2*_ghostLayerWidth" !!!

        assertion(sizeOfPatch(0)==sizeOfPatch(1));

        _postProcessing->mapQuantities(
          offsetOfPatch,
          sizeOfPatch,
          offsetOfPatch + (i-_ghostLayerWidth).convertScalar<double>()* (sizeOfPatch(0)/(numberOfCellsPerAxis)),
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
