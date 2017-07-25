// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www.peano-framework.org
#ifndef _TARCH_LA_LUDECOMPOSITION_H_
#define _TARCH_LA_LUDECOMPOSITION_H_

#include "tarch/la/Vector.h"
#include "tarch/la/Matrix.h"

#include "tarch/Assertions.h"

namespace tarch {
  namespace la {
    /**
     * Performs an in-situ LU-decomposition of the square matrix A. Adds pivot
     * line indices.
     */
    template<int Rows, typename Scalar>
    void lu (
      Matrix<Rows,Rows,Scalar>&  A,
      Vector<Rows,Scalar>&       pivots
    );

    /**
     * In-situ LU without pivoting.
     */
    template<int Rows, typename Scalar>
    void lu (
      Matrix<Rows,Rows,Scalar>&  A
    );

    /**
     * Accepts an upper triangular matrix and a rhs. It then returns the
     * solution x to
     *
     * Rx=f
     *
     * i.e. x=R^{-1}f
     */
    template<int Rows, typename Scalar>
    Vector<Rows,Scalar> backSubstitution(
      const Matrix<Rows,Rows,Scalar>&  R,
      const Vector<Rows,Scalar>&       f
    );
  }
}


#include "tarch/la/LUDecomposition.cpph"


#endif
