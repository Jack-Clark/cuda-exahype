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
 
#ifndef _EXAHYPE_PLOTTERS_FINITE_VOLUMES_2_VTK_ASCII_H_
#define _EXAHYPE_PLOTTERS_FINITE_VOLUMES_2_VTK_ASCII_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

namespace exahype {
  namespace plotters {
    class FiniteVolumes2VTKAscii;
  }
}

class exahype::plotters::FiniteVolumes2VTKAscii: public exahype::plotters::Plotter::Device {
 private:
  int           _fileCounter;
  std::string   _filename;
  int           _numberOfCellsPerAxis;
  int           _ghostLayerWidth;
  int           _solverUnknowns;
  int           _writtenUnknowns;
  std::string   _select;

  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestLeftBottomFront;
  tarch::la::Vector<DIMENSIONS, double>  _regionOfInterestRightTopBack;

  tarch::plotter::griddata::blockstructured::PatchWriterUnstructured*
      _patchWriter;
  tarch::plotter::griddata::blockstructured::PatchWriter::SinglePatchWriter*
      _gridWriter;

  tarch::plotter::griddata::Writer::CellDataWriter*  _timeStampDataWriter;
  tarch::plotter::griddata::Writer::CellDataWriter*  _cellDataWriter;

  static tarch::logging::Log _log;

 public:
  FiniteVolumes2VTKAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing,const int ghostLayerWidth);
  virtual ~FiniteVolumes2VTKAscii();

  virtual void init(const std::string& filename, int numberOfCellsPerAxis, int unknowns, int writtenUnknowns, const std::string& select);

  static std::string getIdentifier();

  void plotPatch(const int cellDescriptionsIndex, const int element) override;

  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);

  virtual void startPlotting( double time );
  virtual void finishPlotting();
};

#endif
