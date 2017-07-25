#include "exahype/plotters/ADERDG2ProbeAscii.h"

#include <fstream>
#include <cstring>

#include "peano/utils/Loop.h"

#include "kernels/DGBasisFunctions.h"

#include "exahype/solvers/ADERDGSolver.h"

tarch::logging::Log exahype::plotters::ADERDG2ProbeAscii::_log( "exahype::plotters::ADERDG2ProbeAscii" );


exahype::plotters::ADERDG2ProbeAscii::ADERDG2ProbeAscii(exahype::plotters::Plotter::UserOnTheFlyPostProcessing* postProcessing):
  Device(postProcessing),
  _out(nullptr) {
}


exahype::plotters::ADERDG2ProbeAscii::~ADERDG2ProbeAscii() {
  if (_out!=nullptr) {
    _out->close();
    _out=nullptr;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::startPlotting( double time ) {
  // In the very first time step, the
  if (time==std::numeric_limits<double>::max() ) {
    _time = 0.0;
  }
  else {
    _time = time;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::finishPlotting() {
  if (_out!=nullptr && *_out) {
    (*_out) << std::endl;
  }
}


void exahype::plotters::ADERDG2ProbeAscii::init(const std::string& filename, int orderPlusOne, int unknowns, int writtenUnknowns, const std::string& select) {
  _order           = orderPlusOne-1;
  _solverUnknowns  = unknowns;
  _writtenUnknowns = writtenUnknowns;
  _select          = select;
  _filename        = filename;
  _time            = 0.0;

  _x(0) = Parser::getValueFromPropertyString( select, "x" );
  _x(1) = Parser::getValueFromPropertyString( select, "y" );
  #if DIM3
  _x(2) = Parser::getValueFromPropertyString( select, "z" );
  #endif

  logDebug( "init(...)", "probe at location " << _x );

  if (!tarch::la::equals(_x,_x)) {
    logError( "init(...)", "Probe location is invalid." );
  }
}


void exahype::plotters::ADERDG2ProbeAscii::openOutputStream() {
  if (_out == nullptr) {
    if (tarch::la::oneEquals(_x,std::numeric_limits<double>::quiet_NaN())) {
      logError( "init(...)", "probe requires valid x, y (and z) coordinates in select statement. No plot written as plot location has been " << _x );
    }
    else {
      _out = new std::ofstream;

      std::ostringstream outputFilename;
      outputFilename << _filename
                   #ifdef Parallel
	           << "-rank-" << tarch::parallel::Node::getInstance().getRank()
                   #endif
		           << ".probe";
      _out->open( outputFilename.str() );
      
      // See issue #47 for discussion whether to quit program on failure
      if(*_out && _out->fail()) {
         logError("openOutputStream(...)", "Could not open file '" << outputFilename.str() << "': " << strerror(errno));
	 exit(-2);
      }

      if (*_out) {
        (*_out) << "# plot-time, real-time";
        for (int unknown=0; unknown < _writtenUnknowns; unknown++) {
          std::ostringstream identifier;
          identifier << "Q" << unknown;
          (*_out) << "," << identifier.str();
        }
        (*_out) << std::endl;
      }
    }
  }

  if (_out!=nullptr && *_out  ) {
    (*_out) << _time;
  }
}


std::string exahype::plotters::ADERDG2ProbeAscii::getIdentifier() {
  return "probe::ascii";
}

void exahype::plotters::ADERDG2ProbeAscii::plotPatch(const int cellDescriptionsIndex, const int element) {
  auto& aderdgCellDescription = exahype::solvers::ADERDGSolver::getCellDescription(cellDescriptionsIndex,element);

  if (aderdgCellDescription.getType()==exahype::solvers::ADERDGSolver::CellDescription::Type::Cell) {
    double* solverSolution = DataHeap::getInstance().getData(aderdgCellDescription.getSolution()).data();

    plotPatch(
        aderdgCellDescription.getOffset(),
        aderdgCellDescription.getSize(), solverSolution,
        aderdgCellDescription.getCorrectorTimeStamp());
  }
}

void exahype::plotters::ADERDG2ProbeAscii::plotPatch(
  const tarch::la::Vector<DIMENSIONS, double>& offsetOfPatch,
  const tarch::la::Vector<DIMENSIONS, double>& sizeOfPatch, double* u,
  double timeStamp
) {
  if (
    tarch::la::allSmaller(offsetOfPatch,_x)
    &&
    tarch::la::allGreater(offsetOfPatch+sizeOfPatch,_x)
  ) {
    // lazy opening
    openOutputStream();


    double* interpoland = new double[_solverUnknowns];
    double* value       = _writtenUnknowns==0 ? nullptr : new double[_writtenUnknowns];

    for (int unknown=0; unknown < _solverUnknowns; unknown++) {
      interpoland[unknown] = kernels::interpolate(
        offsetOfPatch.data(),
        sizeOfPatch.data(),
        _x.data(),
        _solverUnknowns,
        unknown,
        _order,
        u
      );
    }

    _postProcessing->mapQuantities(
      offsetOfPatch,
      sizeOfPatch,
      _x,
      tarch::la::Vector<DIMENSIONS, int>(0),
      interpoland,
      value,
      timeStamp
    );

    if (_out!=nullptr && *_out) {
      (*_out) << ", " << timeStamp;
      for (int i=0; i<_writtenUnknowns; i++) {
        (*_out) << ", " << value[i];
      }
    }

    delete[] interpoland;
    delete[] value;
  }
}
