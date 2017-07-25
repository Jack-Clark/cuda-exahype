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
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_VTK_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_VTK_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/UnstructuredGridWriter.h"

namespace exahype {
  namespace plotters {
    class ADERDG2LegendreVTK;

    class ADERDG2LegendreVerticesVTKAscii;
    class ADERDG2LegendreVerticesVTKBinary;
    class ADERDG2LegendreCellsVTKAscii;
    class ADERDG2LegendreCellsVTKBinary;
  }
}

/**
 * Common VTK class. Usually not used directly but through one of the subclasses.
 */
class exahype::plotters::ADERDG2LegendreVTK: public exahype::plotters::Plotter::Device {
 private:
  int           _fileCounter;
  const bool    _isBinary;
  const bool    _plotCells;
  std::string   _filename;
  int           _order;
  int           _solverUnknowns;
  int           _writtenUnknowns;
  std::string   _select;


  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestLeftBottomFront;
  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestRightTopBack;

  tarch::plotter::griddata::unstructured::UnstructuredGridWriter*                    _gridWriter;

  tarch::plotter::griddata::unstructured::UnstructuredGridWriter::VertexWriter*              _vertexWriter;
  tarch::plotter::griddata::unstructured::UnstructuredGridWriter::CellWriter*                _cellWriter;

  tarch::plotter::griddata::Writer::VertexDataWriter*  _vertexTimeStampDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*    _cellTimeStampDataWriter;
  tarch::plotter::griddata::Writer::VertexDataWriter*  _vertexDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*    _cellDataWriter;

  void writeTimeStampDataToPatch( double timeStamp, int vertexIndex, int cellIndex );

  void plotVertexData(
    int firstVertexIndex,
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp
  );

  void plotCellData(
    int firstCellIndex,
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp
  );

  std::pair<int,int> plotLegendrePatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch
  );
 public:
  ADERDG2LegendreVTK(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool plotCells);
  virtual ~ADERDG2LegendreVTK();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  void plotPatch(const int cellDescriptionsIndex, const int element) override;

  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);

  virtual void startPlotting( double time );
  virtual void finishPlotting();
};


class exahype::plotters::ADERDG2LegendreVerticesVTKAscii: public exahype::plotters::ADERDG2LegendreVTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreVerticesVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreVerticesVTKBinary: public exahype::plotters::ADERDG2LegendreVTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreVerticesVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreCellsVTKAscii: public exahype::plotters::ADERDG2LegendreVTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreCellsVTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


class exahype::plotters::ADERDG2LegendreCellsVTKBinary: public exahype::plotters::ADERDG2LegendreVTK {
  public:
    static std::string getIdentifier();
    ADERDG2LegendreCellsVTKBinary(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
};


#endif
