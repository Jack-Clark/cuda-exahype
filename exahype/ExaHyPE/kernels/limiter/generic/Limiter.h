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

#ifndef _EXAHYPE_KERNELS_LIMITER_GENERIC_H_
#define _EXAHYPE_KERNELS_ADERDG_GENERIC_H_

#include <algorithm>
#include <stdexcept>
#include <stdlib.h>

#include "../../LimiterProjectionMatrices.h"
#include "../../GaussLegendreQuadrature.h"
#include "../../KernelUtils.h"

#include "peano/utils/Globals.h"

namespace kernels {
namespace limiter {
namespace generic {
namespace c {

/**
 * \brief Projection ADERDG -> FV
 *
 * Projects the ADERDG solution onto
 * the finite volumes limiter space.
 *
 * \param[in] basisSize The size of the ADER-DG basis per coordinate axis (order+1).
 * \param[in] ghostLayerWidth The ghost layer width in cells of the finite volumes patch
 */
void projectOnFVLimiterSpace(const double* const luh, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const lim);
/**
 * \brief Projection FV -> ADERDG
 *
 * Projects the finite volumes limiter solution onto
 * the DG space.
 *
 * \param[in] basisSize The size of the ADER-DG basis per coordinate axis (order+1)
 * \param[in] ghostLayerWidth The ghost layer width in cells of the finite volumes patch.
 */
void projectOnDGSpace(const double* const lim, const int numberOfVariables, const int basisSize, const int ghostLayerWidth, double* const luh);

// Get the local min/max from the DG and Gauss Lobatto nodes
void findCellLocalMinAndMax(const double* const luh, const int numberOfVariables, const int basisSize, double* const localMinPerVariables, double* const localMaxPerVariables);

/**
 * Find the minimum and maximum per variable in the limiter solution.
 *
 * We need this function to compute the minimum and maximum
 * values for cells that do not hold a valid ADER-DG
 * solution (troubled cells) and their neighbours.
 * See SUBROUTINE GetMinMax in file
 * ADERDG_Limiter_3D/Limiter.f90.
 *
 * \param[in] basisSize The size of the ADER-DG basis per coordinate axis (order+1)
 * \param[in] ghostLayerWidth The ghost layer width in cells of the finite volumes patch.
 */
void findCellLocalLimiterMinAndMax(const double* const lim, const int numberOfVariables, const int basisSize, const int ghostLayerWidth,
                                   double* const localMinPerVariables, double* const localMaxPerVariables);

/**
 * Similar to ::discreteMaximumPrinciple
 * but writes back the computed cell-local min and max to the
 * boundaryMinPerVariables and boundaryMaxPerVariables arrays.
 */
bool discreteMaximumPrincipleAndMinAndMaxSearch(const double* const luh, const int numberOfVariables, const int basisSize,
                    const double DMPMaximumRelaxationParameter,const double DMPDifferenceScaling,
                    double* boundaryMinPerVariables, double* boundaryMaxPerVariables);

/**
 * Returns true if the nodal solution degrees of freedom
 * satisfy a discrete maximum principle.
 *
 * \note[24/11/16]
 * We currently abuse the term Voronoi neighbour for direct neighbour.
 *
 * \param[in] luh                           The nodal solution degrees of freedom
 * \param[in] numberOfVariables             The number of variables of the solved PDE
 * \param[in] basisSize                     The basis size of the ADER-DG discretisation (approx. order + 1).
 * \param[in] DMPMaximumRelaxationParameter The relaxation parameter for the discrete maximum principle (DMP).
 * \param[in] DMPDifferenceScaling          The difference scaling factor for the discrete maximum principle (DMP).
 * \param[in] boundaryMinPerVariables       An array of size \p numberOfVariables times DIMENSIONS_TIMES_TWO
 *                                          containing the minimum values per variable of the current cell
 *                                          and its neighbour at the particular face. Together these values
 *                                          can be used to compute the Voronoi maximum per variable.
 * \param[in] boundaryMinPerVariables       An array of size \p numberOfVariables times DIMENSIONS_TIMES_TWO
 *                                          containing the minimum values per variable of the current cell
 *                                          and its neighbour at the particular face. Together these values
 *                                          can be used to compute the Voronoi maximum per variable.
 *
 * <h2>Background</h2>
 * A candidate solution \f$ u^{*}_h(x,t^{n+1}) \f$ is said to satisfy
 * the discrete maximum principle if it satisfies a relaxed
 * maximum principle of the form
 *
 * \f[
 *   \min_{y \in V_i} (u_h(y,t^n)) - \delta \leq \, u^{*}_h(x,t^{n+1}) \leq \, \max_{y \in V_i} (u_h(y,t^n)) + \delta,
 *   \;\forall \x \in T_i
 * \f]
 *
 * for every element \f$ T_i \f$ in the mesh. Above, $\f$ V_i \f$ is a set containing the Voronoi neighbours
 * of element \f$ T_i \f$ and the \f$ T_i \f$ itself.
 * The relaxation parameter \f$ \delta \f$ is computed according to:
 *
 * \f[
 *  \delta = \max \left( \delta_0,\, \epsilon \cdot \left( \max_{y \in V_i} (u_h(y,t^n)) - \min_{y \in V_i} (u_h(y,t^n)) \right) \right),
 * \f]
 * with \f$ \delta_0 \f$ denoting the maximum relaxation parameter we want to allow,
 * and \epsilon scales the difference of Voronoi maximum and minimum.
 *
 * See doi:10.1016/j.jcp.2014.08.009 for more details.
 */
bool discreteMaximumPrinciple(const double* const luh, const int numberOfVariables, const int basisSize,
                    const double DMPMaximumRelaxationParameter,const double DMPDifferenceScaling,
                    const double* const boundaryMinPerVariables, const double* const boundaryMaxnPerVariables);

// TODO(Dominic): @JM: We have to do a rollback in every neighbour cell of the troubled cells. Furthermore, the
// troubled cells are not that many compared to the non-troubled ones. Thus, I decided to get
// rid of the solution anticipation. I am sorry for the confusion.
//bool isTroubledCell(const double* const luh, const double* const lduh, const double dt, const int numberOfVariables, const int basisSize, const double* const boundaryMinPerVariables, const double* const boundaryMaxPerVariables);

//************************
//*** Helper functions ***
//************************

//inline double anticipateLuh(const double* const luh, const double* const lduh, const double dt, const int order, const int idx, const int x, const int y, const int z) {
//  double weight =
//  #if DIMENSIONS == 3
//  kernels::gaussLegendreWeights[order][z] *
//  #endif
//  kernels::gaussLegendreWeights[order][y] * kernels::gaussLegendreWeights[order][x];
//
//  return (luh[idx] + dt / weight * lduh[idx]); // TODO(Dominic): The compiler might not able to optimise for the dt=0 case.
//}

inline int getBasisSizeLim(const int basisSize) {
  return 2*(basisSize-1)+1;
}

//*************************
//*** Private functions ***
//*************************

//Projection ADERDG -> Gauss-Lobatto, for test only
double* getGaussLobattoData(const double* const luh, const int numberOfVariables, const int basisSize); 

void compareWithADERDGSolutionAtGaussLobattoNodes(const double* const luh, const int numberOfVariables, const int basisSize, double* const min, double* const max);

} // namespace c
} // namespace generic
} // namespace limiter
} // namespace kernel


#endif //_EXAHYPE_KERNELS_LIMITER_GENERIC_H_
