/** \file IOperator.hpp
 *  \brief Interface for a general spectral operator implementation (direct and quasi inverses)
 */

#ifndef IOPERATOR_HPP
#define IOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Simulation/Enums/BoundaryConditions.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Interface for a general spectral operator implementation (direct and quasi inverses)
    */
   class IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param idx        Index of assigned parameter
          * @param specSize   Size of the spectral space
          * @param specIdx    Spectral indexes
          */
         SpectralOperatorBase(const int idx, const int specSize, const ArrayI& specIdx);

         /**
          * @brief Empty Destructor
          */
         virtual ~SpectralOperatorBase();

         /**
          * @brief Reset spectral operators and set new dimension and spectral indexes
          *
          * @param specSize   New size of the spectral space
          * @param specIdx    New set of spectral indexes
          */
         void reset(const int specSize, const ArrayI& specIdx);

         /**
          * @brief Get the (qasi) identity matrix of order p
          *
          * The 0th order quasi identity is the standard identity matrix. For positiv p the top rows are removed
          * and for negative p the bottom rows are removed
          *
          * @param p Order of the quasi identity
          */
         virtual SparseMatrix id(const int p = 0);

         /**
          * @brief Get the derivative of order p
          *
          * @param nBC  Number of boundary
          * @param p    Order of the derivative
          */
         virtual SparseMatrix diff(const int nBC, const int p) = 0;

         /**
          * @brief Get the quasi inverse derivative for D^-p D^q
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         virtual SparseMatrix qDiff(const int p, const int q) = 0;
         
      protected:
         /**
          * @brief Get size of spectral space
          */
         int polyN() const;

         /**
          * @brief Convert tau lines into sparse matrix
          *
          * @param lines   Tau lines
          * @param atTop   Put tau lines at top? (default: true)
          */
         DecoupledZSparse  createSparseTau(const DecoupledZMatrix& lines, const bool atTop = true) const;

      private:
         /**
          * @brief Size of the polynomial basis
          */
         int mPolyN;
   };

}
}

#endif // IOPERATOR_HPP
