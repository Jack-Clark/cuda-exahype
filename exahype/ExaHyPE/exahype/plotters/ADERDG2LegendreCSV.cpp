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
 
#include "ADERDG2LegendreCSV.h"
#include "tarch/parallel/Node.h"

#include "kernels/DGMatrices.h"
#include "kernels/GaussLegendreQuadrature.h"
#include "kernels/DGBasisFunctions.h"

#include "peano/utils/Loop.h"

#include "exahype/solvers/ADERDGSolver.h"

#include <iomanip>

#define eps 1e-13
using namespace std;



std::string exahype::plotters::ADERDG2LegendreCSV::getIdentifier() {
  return "csv::Legendre::nodes::ascii";
}

exahype::plotters::ADERDG2LegendreCSV::ADERDG2LegendreCSV(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
    ADERDG2LegendreCSV(postProcessing,false,false) {
}

exahype::plotters::ADERDG2LegendreCSV::ADERDG2LegendreCSV(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing, bool isBinary, bool plotCells):
  Device(postProcessing),
  _fileCounter(-1),
  _isBinary(isBinary),
  _plotCells(plotCells) {
    
    
}


void exahype::plotters::ADERDG2LegendreCSV::init(
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
  
#ifndef Dim2
    cout << "ADERDG2LegendreCSV::init(...)" << "Plotter is not yet implemented for 3d" << endl;
    exit(-42);
  
#endif  
}


void exahype::plotters::ADERDG2LegendreCSV::startPlotting( double time ) {
  _fileCounter++;


  std::ostringstream fileName;
  fileName << _filename << "-" << _fileCounter << ".csv";  
  ofs.open(fileName.str(), std::ofstream::out);
  
  //Header
  // ofs << "x, y, timestamp";
  // for (int unknown=0; unknown < _solverUnknowns; unknown++) {
    // ofs << ", Q" << unknown;
  // }
  // ofs << endl;
  

  if (_writtenUnknowns>0) {

  }

  _postProcessing->startPlotting( time );
}


void exahype::plotters::ADERDG2LegendreCSV::finishPlotting() {
  

  ofs.close();

  _postProcessing->finishPlotting();

  if (_writtenUnknowns>0) {

  }
}


exahype::plotters::ADERDG2LegendreCSV::~ADERDG2LegendreCSV() {
}


void exahype::plotters::ADERDG2LegendreCSV::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2LegendreCSV::plotPatch(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    double* u,
    double timeStamp) {
      
      
 
  // VV-Todo assert(dim==2)
  // cout << offsetOfPatch(0) << endl;

  // int treeDepth = log(sizeOfPatch(0))/log(1./3.) + eps;
  // double increment = pow(1./3., treeDepth);
  // int xIndex = offsetOfPatch(0)/increment + eps;
  // int yIndex = offsetOfPatch(1)/increment + eps;
  // int elementsPerAxis = pow(3., treeDepth) + eps;

  // cout << offsetOfPatch << endl;
  // cout << sizeOfPatch << endl;
  // cout << "treeDepth= " << treeDepth << endl;
  // cout << "elementsPerAxis= " << elementsPerAxis << endl;
  // cout << "increment= " << increment << endl;
  // cout << "yIndex= " << yIndex << endl;
  // cout << "xIndex= " << xIndex << endl;
  // cout << "index= " << yIndex*elementsPerAxis + xIndex << endl;

  
  for (int x=0; x < _order + 1; x++){
    for (int y=0; y < _order + 1; y++){
//      double weight = kernels::gaussLegendreWeights[_order][x] *
//                      kernels::gaussLegendreWeights[_order][y] *
//                      sizeOfPatch(0)*sizeOfPatch(1); // TODO(Dominic): Unused
                      
      
      //Data
      double dx = offsetOfPatch(0) + sizeOfPatch(0) * kernels::gaussLegendreNodes[_order][x];
      double dy = offsetOfPatch(1) + sizeOfPatch(1) * kernels::gaussLegendreNodes[_order][y];
      
      ofs << std::setprecision(14);
      ofs << dx << ", " << dy << ", " << timeStamp;
      
      for (int unknown=0; unknown < _solverUnknowns; unknown++) {
        int idx = _solverUnknowns * ((_order + 1)*y + x) + unknown;
        
        // cout << "u[" << idx << "]= " << u[idx] << endl;
        ofs << ", " << u[idx];
      }
      ofs << endl;
    }
  }
  // cout << endl;

}
