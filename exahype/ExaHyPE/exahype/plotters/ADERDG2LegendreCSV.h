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
 
#ifndef _EXAHYPE_PLOTTERS_ADERDG_2_CSV_H_
#define _EXAHYPE_PLOTTERS_ADERDG_2_CSV_H_

#include "exahype/plotters/Plotter.h"

#include "tarch/plotter/griddata/blockstructured/PatchWriterUnstructured.h"
#include "tarch/plotter/griddata/unstructured/UnstructuredGridWriter.h"

#include <fstream> 

namespace exahype {
  namespace plotters {
    class ADERDG2LegendreCSV;

  }
}

/**
 * todo
 */
class exahype::plotters::ADERDG2LegendreCSV: public exahype::plotters::Plotter::Device {
 private:
  int           _fileCounter;
  const bool    _isBinary;
  const bool    _plotCells;
  std::string   _filename;
  int           _order;
  int           _solverUnknowns;
  int           _writtenUnknowns;
  std::string   _select;


  std::ofstream ofs;

 public:
  static std::string getIdentifier();
  ADERDG2LegendreCSV(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing);
  ADERDG2LegendreCSV(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool plotCells);
  virtual ~ADERDG2LegendreCSV();

  virtual void init(const std::string& filename, int orderPlusOne, int solverUnknowns, int writtenUnknowns, const std::string& select);

  void plotPatch(const int cellDescriptionsIndex, const int element) override;

  void plotPatch(
      const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
      const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
      double timeStamp);

  virtual void startPlotting( double time );
  virtual void finishPlotting();
};
#endif
