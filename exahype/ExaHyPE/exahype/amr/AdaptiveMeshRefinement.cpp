#include "exahype/amr/AdaptiveMeshRefinement.h"

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

tarch::la::Vector<DIMENSIONS,int> exahype::amr::computeSubcellIndex(
      const tarch::la::Vector<DIMENSIONS,double>& childOffset,
      const tarch::la::Vector<DIMENSIONS,double>& childSize,
      const tarch::la::Vector<DIMENSIONS,double>& parentOffset) {
    tarch::la::Vector<DIMENSIONS,int> subcellIndex;
    for (int xi = 0; xi < DIMENSIONS; ++xi) {
      assertion((childOffset(xi) - parentOffset(xi)) >= 0);
      subcellIndex[xi] = tarch::la::round(
          (childOffset(xi) - parentOffset(xi))/childSize(xi));
    }

    return subcellIndex;
  }

tarch::la::Vector<DIMENSIONS-1,int> exahype::amr::getSubfaceIndex(
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int d) {
  tarch::la::Vector<DIMENSIONS-1,int> subfaceIndex;

  int i = 0;
  for (int j = 0; j < DIMENSIONS; j++) {
    if (j != d) {
      subfaceIndex[i] = subcellIndex[j];
      i++;
    }
  }

  return subfaceIndex;
}

bool exahype::amr::onBoundaryOfParent(
    const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
    const int levelDelta){
  for (int d = 0; d < DIMENSIONS; d++) {
    // Check if cell is at "left" or "right" d face of parent
    if (subcellIndex[d] == 0 ||
        subcellIndex[d] == tarch::la::aPowI(levelDelta,3)-1) {
      return true;
    }
  }
  return false;
}
