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

#ifndef _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_H_
#define _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"

namespace exahype {
  namespace plotters {
    class LimitingADERDG2CartesianVTK;

    class LimitingADERDG2CartesianVerticesVTKAscii;
    class LimitingADERDG2CartesianVerticesVTKBinary;
    class LimitingADERDG2CartesianCellsVTKAscii;
    class LimitingADERDG2CartesianCellsVTKBinary;
  }
}

/**
 * Common VTK class. Usually not used directly but through one of the subclasses.
 */
class exahype::plotters::LimitingADERDG2CartesianVTK: public exahype::plotters::Plotter::Device {
private:
  int           _fileCounter;
  const bool    _isBinary;
  const bool    _plotCells;
  std::string   _filename;
  int           _order;
  int           _solverUnknowns;
  int           _writtenUnknowns;
  /**
   * The ghost layer width the finite volumes patch is using.
   */
  const int     _ghostLayerWidth;
  std::string   _select;

  static tarch::logging::Log _log;

  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestLeftBottomFront;
  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestRightTopBack;

  tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter* _gridWriter;
  tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*        _patchWriter;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _vertexDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _cellDataWriter;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _timeStampVertexDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _timeStampCellDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*                          _cellLimiterStatusWriter;
  tarch::plotter::griddata::Writer::VertexDataWriter*                        _vertexLimiterStatusWriter;

  void writeTimeStampDataToADERDGPatch( double timeStamp, int vertexIndex );

  void plotVertexData(
      int firstVertexIndex,
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      double* u,
      double timeStamp,
      const int limiterStatus
  );

  void plotCellData(
      int firstCellIndex,
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
      double* u,
      double timeStamp,
      const int limiterStatusAsInt
  );
public:
  LimitingADERDG2CartesianVTK(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                              const int ghostLayerWidth,const bool isBinary, const bool plotCells);
  virtual ~LimitingADERDG2CartesianVTK();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  void plotPatch(const int cellDescriptionsIndex, const int element) override;

  /**
   * Plot an ADER-DG solution.
   */
  void plotADERDGPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp,
      const int limiterStatusAsInt);

  /**
   * Plot a finite volumes solution.
   */
  void plotFiniteVolumesPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);

  void startPlotting( double time ) override;
  void finishPlotting() override;
};

class exahype::plotters::LimitingADERDG2CartesianVerticesVTKAscii: public exahype::plotters::LimitingADERDG2CartesianVTK {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                             const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianVerticesVTKBinary: public exahype::plotters::LimitingADERDG2CartesianVTK {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                              const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianCellsVTKAscii: public exahype::plotters::LimitingADERDG2CartesianVTK {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                          const int ghostLayerWidth);
};


class exahype::plotters::LimitingADERDG2CartesianCellsVTKBinary: public exahype::plotters::LimitingADERDG2CartesianVTK {
  public:
    static std::string getIdentifier();
    LimitingADERDG2CartesianCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,
                                           const int ghostLayerWidth);
};

#endif // _EXAHYPE_PLOTTERS_LIMITING_ADERDG_2_CARTESIAN_VTK_H_
