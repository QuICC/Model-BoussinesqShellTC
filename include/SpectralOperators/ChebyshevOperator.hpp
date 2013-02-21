/** \file ChebyshevOperator.hpp
 *  \brief Implementation of the spectral operators for the chebyshev basis
 */

#ifndef CHEBYSHEVOPERATOR_HPP
#define CHEBYSHEVOPERATOR_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/General/Typedefs.hpp"
#include "SpectralOperators/IOperatorBase.hpp"
#include "Simulation/Enums/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * \brief Implementation of the spectral operators for the chebyshev basis
    */
   class ChebyshevOperator: public IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param idx     Index of assigned parameter
          * @param polyN   Size of the polynomial basis
          * @param specIdx Spectral indexes
          */
         ChebyshevOperator(const int idx, const int polyN, const ArrayI& specIdx);

         /**
          * @brief Empty Destructor
          */
         ~ChebyshevOperator();

         /**
          * @brief Get the derivative operator of order p
          *
          * @param nBC  Number of boundary
          * @param p    Order of the derivative
          */
         SparseMatrix diff(const int nBC, const int p);

         /**
          * @brief Get the quasi inverse derivative for D^-p D^q
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         SparseMatrix qDiff(const int p, const int q);

         /**
          * @brief Convert tau lines to complete sparse matrix
          *
          * @param bcId    List of BC to implement
          * @param atTop   Set tau lines at top? (default: true)
          */
         DecoupledZSparse tau(const std::map<BoundaryConditions::Id,BoundaryConditions::Position>& bcId, const bool atTop = true);
         
      protected:

      private:
         /**
          * @brief C factor
          *
          * @param n Order of polynial
          */
         EPMFloat c(const int n) const;

         /**
          * @brief Build the derivative matrix
          *
          * @param mat Storage to build the matrix into
          */
         void buildDerivative(SparseMatrix& mat) const;

         /**
          * @brief Build the inverse matrix
          *
          * @param mat Storage to build the matrix into
          */
         void buildInverse(SparseMatrix& mat) const;
   };

}
}

#endif // CHEBYSHEVOPERATOR_HPP
