/** \file UnitOperator.cpp
 *  \brief Source of the implementation of the unit spetral operator and quasi inverses
 */

// System includes
//
#include <assert.h>

// External includes
//

// Class include
//
#include "SpectralOperators/UnitOperator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Spectral {

   UnitOperator::UnitOperator(const int basisN)
      : IOperator(basisN)
   {
   }

   UnitOperator::~UnitOperator()
   {
   }

   SparseMatrix UnitOperator::id(const int p)
   {
      // Create the 1x1 identity
      SparseMatrix idMat(1, 1);

      // Create left quasi identity
      idMat.reserve(idMat.cols());
      for(int j = 0; j < idMat.cols(); ++j)
      {
         // Create column j
         idMat.startVec(j);

         idMat.insertBack(j,j) = 1.0;
      }
      idMat.finalize(); 

      return idMat;
   }

   SparseMatrix UnitOperator::diff(const int nBC, const int p)
   {
      // Check that requested order is possible
      assert(p > 0);

      // 1x1 identity
      return this->id();
   }

   SparseMatrix UnitOperator::qDiff(const int p, const int q)
   {
      // Check that requested order is possible
      assert(p > 0);
      assert(q > 0);
      assert(p > q);

      // 1x1 identity
      return this->id();
   }

}
}
