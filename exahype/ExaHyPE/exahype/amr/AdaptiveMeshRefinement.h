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

#ifndef ADAPTIVEMESHREFINEMENT_H_
#define ADAPTIVEMESHREFINEMENT_H_

#include "peano/utils/Globals.h"

#include "exahype/Cell.h"

#include "multiscalelinkedcell/HangingVertexBookkeeper.h"

namespace exahype {

namespace amr {
  /**
   * Per coordinate direction xi, count the number of shifts
   * of step size \p childSize(xi) necessary to
   * reach \p childOffset from \p parentOffset.
   *
   * \param[in] childOffset  Offset of a child cell.
   * \param[in] childSize    Size of the child cell.
   * \param[in] parentOffset Offset of the parent cell.
   *
   * \see getSubfaceIndex
   */
  tarch::la::Vector<DIMENSIONS,int> computeSubcellIndex(
        const tarch::la::Vector<DIMENSIONS,double>& childOffset,
        const tarch::la::Vector<DIMENSIONS,double>& childSize,
        const tarch::la::Vector<DIMENSIONS,double>& parentOffset);

  /**
   * Collect all the element with index!=d
   * from \p subcellIndex.
   *
   * \see getSubcellIndex
   */
  tarch::la::Vector<DIMENSIONS-1,int> getSubfaceIndex(
        const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
        const int d);

  /**
   * Per coordinate direction xi, check if
   * subcellIndex[xi] is either 0 or 3^levelDelta - 1.
   * If this is the case for at least one subcellIndex[xi]
   * return true. Otherwise return false.
   */
  bool onBoundaryOfParent(
      const tarch::la::Vector<DIMENSIONS, int>& subcellIndex,
      const int levelDelta);

  /**
   * Returns  AugmentationControl::::NextToCell if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Cell,
   * exahype::solvers::Solver::NextToAncestor if the cell has a neighbour of
   * type exahype::records::ADERDGCellDescription:Ancestor or
   * exahype::records::ADERDGCellDescription::EmptyAncestor.
   * Is both the case, this function returns
   * AugmentationControl::::NextToCellOrAncestor.
   * If none of the previous is the case, this function returns
   * AugmentationControl::Default.
   */
  template<class CellDescription,class CellDescriptionsHeap>
  exahype::solvers::Solver::AugmentationControl
  augmentationCriterion(
      const int solverNumber,
      const typename CellDescription::Type type,
      const int level,
      const tarch::la::Vector<THREE_POWER_D, int>&
      neighbourCellDescriptionsIndices) {
    // left,right,front,back,(bottom,top)
#if DIMENSIONS == 2
    constexpr int neighbourPositions[4] = {3, 5, 1, 7};
#else
    constexpr int neighbourPositions[6] = {12, 14, 10, 16, 4, 22};
#endif
    bool nextToAncestor = false;
    bool nextToCell = false;

    for (int i = 0; i < DIMENSIONS_TIMES_TWO; i++) {
      const int neighbourCellDescriptionIndex =
          neighbourCellDescriptionsIndices[neighbourPositions[i]];
      if (CellDescriptionsHeap::getInstance().isValidIndex(neighbourCellDescriptionIndex)) {
        for (CellDescription& pNeighbour : CellDescriptionsHeap::getInstance().getData(
            neighbourCellDescriptionIndex)) {
          if (pNeighbour.getSolverNumber() == solverNumber &&
              pNeighbour.getLevel() == level) {
            switch (pNeighbour.getType()) {
              case CellDescription::Ancestor:
              case CellDescription::EmptyAncestor:
                nextToAncestor = true;
                break;
              case CellDescription::Cell:
                nextToCell = true;
                break;
              default:
                break;
            }
          }
        }
      }
    }

    // NOTE: The order below is important.
    if (nextToCell && nextToAncestor) {
      return exahype::solvers::Solver::AugmentationControl::NextToCellAndAncestor;
    }
    if (nextToAncestor) {
      return exahype::solvers::Solver::AugmentationControl::NextToAncestor;
    }
    if (nextToCell) {
      return exahype::solvers::Solver::AugmentationControl::NextToCell;
    }
    // Erase otherwise.
    return exahype::solvers::Solver::AugmentationControl::Default;
  }

  /*
   * Change the erasing request to a change to descendant request if the coarse grid Cell
   * has children (of type Descendant).
   * Rationale: We cannot directly erase a Cell that has children (of type Descendant).
   *
   * Further, reset the deaugmenting request if a coarse grid Descendant has children
   * (of type Descendant). Rationale:
   * We cannot erase a coarse grid cell that has children (of type Descendant)
   * before erasing the children.
   *
   * \note This method should be called when we enter a fine grid cell,
   * coarseGridCellDescription is a cell description associated with
   * the coarse grid cell.
   *
   * \note A more sophisticated procedure has to performed for the refinement event
   * AugmentationRequested. We need to use the taversal's descend event to handle
   * this event. We thus do not rely on fineGridCell.isRefined() in the previous enterCell event
   * to check if we need to reset the deaugmenting request.
   */
  template<class CellDescription>
  void resetErasingOrDeaugmentingRequestIfParent(
      CellDescription& coarseGridCellDescription) {
    switch (coarseGridCellDescription.getRefinementEvent()) {
      case CellDescription::DeaugmentingRequested:
        assertion1(coarseGridCellDescription.getType()==CellDescription::EmptyDescendant ||
                   coarseGridCellDescription.getType()==CellDescription::Descendant,coarseGridCellDescription.toString());
        coarseGridCellDescription.setRefinementEvent(CellDescription::None);
        break;
      case CellDescription::ErasingRequested:
        assertion1(coarseGridCellDescription.getType()==CellDescription::Cell,
                   coarseGridCellDescription.toString());
        coarseGridCellDescription.setRefinementEvent(CellDescription::ChangeToDescendantRequested);
        break;
      default:
        break;
    }
  }

  /**
   * Determine the position of a Cell or Ancestor with respect
   * to a parent of type Ancestor.
   * The return values subcellPosition.parentCellDescriptionsIndex
   * and subcellPosition.parentElement only
   * hold valid indices >= 0 if we have found a parent of type Ancestor
   * and the cell description itself is of type Cell or Ancestor or EmptyAncestor.
   *
   * Otherwise, subcellPosition.parentIndex holds the value
   * multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
   * subcellPosition.parentElement holds the value exahype::solver::Solvers::NotFound,
   * and subcellPosition.subcellIndex holds undefined values.
   * This applies to the case where the parent index of the
   * Cell or Ancestor is not a valid index and further to
   * the case where we have found a parent of type EmptyAncestor
   * as top most parent.
   *
   * This method is required for preparing cell description types
   * before sending the cell description away to a new worker.
   */
  template <class CellDescription,class CellDescriptionHeap>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfCellOrAncestorOrEmptyAncestor(
      const CellDescription& pChild) {
    // 1. Initialisation.
    exahype::solvers::Solver::SubcellPosition subcellPosition;
    if ((pChild.getType()==CellDescription::Cell
        || pChild.getType()==CellDescription::Ancestor
        || pChild.getType()==CellDescription::EmptyAncestor)
        && CellDescriptionHeap::getInstance().isValidIndex(pChild.getParentIndex())) {
      CellDescription* pParent = 0;
      int parentElement=0;
      for (auto& p : CellDescriptionHeap::getInstance().getData(pChild.getParentIndex())) {
        if (p.getSolverNumber()==pChild.getSolverNumber()) {
          subcellPosition.parentCellDescriptionsIndex = pChild.getParentIndex();
          subcellPosition.parentElement               = parentElement;
          pParent = &p;
          break;
        }
        ++parentElement;
      }

      // 2. If the current parent is an EmptyAncestor,
      // try to determine iteratively the next parent that holds data.
      if (subcellPosition.parentElement!=exahype::solvers::Solver::NotFound) {
        while (pParent->getType()==CellDescription::EmptyAncestor &&
            CellDescriptionHeap::getInstance().isValidIndex(pParent->getParentIndex())) {
          const int currentParentIndex = pParent->getParentIndex();
          int parentElement=0;
          for (auto& p : CellDescriptionHeap::getInstance().getData(currentParentIndex)) {  // Loop over cell descriptions
            if (p.getSolverNumber() == pChild.getSolverNumber()) {
              subcellPosition.parentCellDescriptionsIndex = pParent->getParentIndex();
              subcellPosition.parentElement               = parentElement;
              pParent = &p;
              break;
            }
            ++parentElement;
          }
        }
        assertion2(pParent->getType()==CellDescription::Ancestor ||
                   pParent->getType()==CellDescription::EmptyAncestor,
                   pChild.toString(),pParent->toString());

        if (pParent->getType()==CellDescription::EmptyAncestor) {
          exahype::solvers::Solver::SubcellPosition subcellPositionWithInvalidIndices;
          return subcellPositionWithInvalidIndices;
        }

        subcellPosition.subcellIndex =
            computeSubcellIndex(
                pChild.getOffset(),pChild.getSize(),
                pParent->getOffset());
        subcellPosition.levelDifference =
            pChild.getLevel() - pParent->getLevel();
      }
    }

    return subcellPosition;
  }

  /**
   * Determine the position of a Cell or Ancestor with respect
   * to a parent of type Ancestor.
   * The return values subcellPosition.parentCellDescriptionsIndex
   * and subcellPosition.parentElement only
   * hold valid indices >= 0 if we have found a parent of type Ancestor
   * and the cell description itself is of type Cell or Ancestor.
   *
   * Otherwise, subcellPosition.parentIndex holds the value
   * multiscalelinkedcell::HangingVertexBookkeeper::InvalidAdjacencyIndex,
   * subcellPosition.parentElement holds the value exahype::solver::Solvers::NotFound,
   * and subcellPosition.subcellIndex holds undefined values.
   * This applies to the case where the parent index of the
   * Cell or Ancestor is not a valid index and further to
   * the case where we have found a parent of type EmptyAncestor
   * as top most parent.
   *
   * This method is required for the face data restriction, the
   * volume data restriction, and the FV volume data restriction.
   */
  template <class CellDescription,class CellDescriptionHeap>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfCellOrAncestor(
      const CellDescription& pChild) {
    exahype::solvers::Solver::SubcellPosition subcellPosition;
    if (pChild.getType()==CellDescription::Cell
        || pChild.getType()==CellDescription::Ancestor) {
      subcellPosition = computeSubcellPositionOfCellOrAncestorOrEmptyAncestor
          <CellDescription,CellDescriptionHeap>(pChild);
    }

    return subcellPosition;
  }

  /**
   * Determine the position of a Descendant with respect
   * to a  Cell or Descendant that contains data, i.e.,
   * has at least one neighbour that is a real cell.
   *
   * \note This function only makes sense if the
   * Descendant has a valid parentIndex attribute.
   *
   * \note This method is less complicated that corresponding
   * method for computing the parent of
   * a Cell or Ancestor since a Descendant always(!) has
   * a parent of type Cell or Descendant.
   *
   * This method is required for the face data prolongation, the
   * volume data prolongation, and the FV volume data prolongation.
   */
  template <class CellDescription,class CellDescriptionHeap>
  exahype::solvers::Solver::SubcellPosition
  computeSubcellPositionOfDescendant(
      const CellDescription& pChild) {
    assertion1(pChild.getType()==CellDescription::Descendant,pChild.toString());
    assertion1(CellDescriptionHeap::getInstance().isValidIndex(
        pChild.getParentIndex()),pChild.getParentIndex());

    // 1. Initialisation.
    exahype::solvers::Solver::SubcellPosition subcellPosition;
    subcellPosition.parentCellDescriptionsIndex = pChild.getParentIndex();
    subcellPosition.parentElement=exahype::solvers::Solver::NotFound;
    CellDescription* pParent = 0;
    int parentElement=0;
    for (auto& p : CellDescriptionHeap::getInstance().getData(
        pChild.getParentIndex())) {
      if (p.getSolverNumber()==pChild.getSolverNumber()) {
        subcellPosition.parentElement = parentElement;
        pParent = &p;
        break;
      }
      ++parentElement;
    }
    // Descendant pChild must always have
    // a parent in the parent's cell description array.
    assertion1(pParent!=0,pChild.toString());

    // 2. If the current parent is an EmptyAncestor,
    // try to determine iteratively the next parent that holds data.
    while (pParent->getType()==CellDescription::EmptyDescendant) {
      const int currentParentIndex =
          pParent->getParentIndex();
      assertion1(CellDescriptionHeap::getInstance().isValidIndex(
          currentParentIndex),currentParentIndex); // Must always hold if the current parent is an (Empty)Descendant
      int parentElement=0;
      for (auto& p : CellDescriptionHeap::getInstance().getData(currentParentIndex)) {
        if (p.getSolverNumber()==pChild.getSolverNumber()) {
          subcellPosition.parentCellDescriptionsIndex = pParent->getParentIndex();
          subcellPosition.parentElement               = parentElement;
          pParent = &p;
          break;
        }
        ++parentElement;
      }
    }
    assertion(pParent->getType() == CellDescription::Descendant ||
              pParent->getType() == CellDescription::Cell);

    subcellPosition.subcellIndex =
            computeSubcellIndex(
                pChild.getOffset(),pChild.getSize(),
                pParent->getOffset());
    subcellPosition.levelDifference =
        pChild.getLevel() - pParent->getLevel();

    return subcellPosition;
  }

}  // namespace amr
}  // namespace exahype



#endif /* ADAPTIVEMESHREFINEMENT_H_ */
