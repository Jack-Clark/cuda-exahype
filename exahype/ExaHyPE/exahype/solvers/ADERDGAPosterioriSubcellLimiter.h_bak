/*
 * ADERDGAPosterioriSubcellLimiter.h
 *
 *  Created on: 6 Oct 2016
 *      Author: dominic
 */

#ifndef ADERDGAPOSTERIORISUBCELLLIMITER_H_
#define ADERDGAPOSTERIORISUBCELLLIMITER_H_

#include "CellWiseCoupling.h"

#include "tarch/logging/Log.h"

namespace exahype {
namespace solvers {

class ADERDGAPosterioriSubcellLimiter;

} /* namespace solvers */
} /* namespace exahype */

class exahype::solvers::ADERDGAPosterioriSubcellLimiter : public exahype::solvers::CellWiseCoupling {
private:
  static tarch::logging::Log _log;

  /**
   * Element index of the ADER-DG solver whose solution is to limit in
   * registry exahype::solvers::RegisteredSolvers.
   */
  const int  _aderdgSolverNumber;
  /**
   * Element index of the Finite Volumes solver used
   * for the subcell limiting in registry
   * exahype::solvers::RegisteredSolvers.
   */
  const int  _finiteVolumesSolverNumber;

public:
  ADERDGAPosterioriSubcellLimiter(int aderdgSolverNumber,int finiteVolumesSolverNumber);
  virtual ~ADERDGAPosterioriSubcellLimiter() {};

  // Disallow copy and assignment
  ADERDGAPosterioriSubcellLimiter(const ADERDGAPosterioriSubcellLimiter& other) = delete;
  ADERDGAPosterioriSubcellLimiter& operator=(const ADERDGAPosterioriSubcellLimiter& other) = delete;

  /**
   * TODO(Dominic): Add docu.
   */
  void coupleFirstTime(
      const int cellDescriptionsIndex,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;

  /**
   * TODO(Dominic): Add docu.
   */
  void couple(
      const int cellDescriptionsIndex,
      exahype::Vertex* const fineGridVertices,
      const peano::grid::VertexEnumerator& fineGridVerticesEnumerator) override;
};

#endif /* ADERDGAPOSTERIORISUBCELLLIMITER_H_ */
