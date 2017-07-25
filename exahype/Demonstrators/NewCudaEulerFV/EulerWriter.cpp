#include "EulerWriter.h"


NewCudaEulerFV::EulerWriter::EulerWriter(MyCudaEulerSolver&  solver) {
  // @todo Please insert your code here
}


NewCudaEulerFV::EulerWriter::~EulerWriter() {
  // @todo Please insert your code here
}


void NewCudaEulerFV::EulerWriter::startPlotting(double time) {
  // @todo Please insert your code here
}


void NewCudaEulerFV::EulerWriter::finishPlotting() {
  // @todo Please insert your code here
}


void NewCudaEulerFV::EulerWriter::mapQuantities(
    const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch,
    const tarch::la::Vector<DIMENSIONS, double>& x,
    const tarch::la::Vector<DIMENSIONS, int>&    pos,
    double* Q,
    double* outputQuantities,
    double timeStamp
) {
  for (int i=0; i<5; i++){ 
    outputQuantities[i] = Q[i];
  }
}


