/** 
 * @file IOperator.hpp
 * @brief Interface for a general spectral operator implementation (direct and quasi inverses)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
          * @param basisN   Size of the spectral basis
          */
         IOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         virtual ~IOperator();

         /**
          * @brief Reset spectral operators and set new dimension
          *
          * @param basisN   New size of the spectral basis
          */
         void reset(const int basisN);

         /**
          * @brief Get the (quasi) identity matrix of order p
          *
          * The 0th order quasi identity is the standard identity matrix. For positiv p the top rows are removed
          * and for negative p the bottom rows are removed
          *
          * @param p Order of the quasi identity
          */
         virtual SparseMatrix id(const int p = 0) const;

         /**
          * @brief Get the shifted (quasi) identity matrix of order p
          *
          * The 0th order shifted quasi identity is the standard identity matrix. For positiv p the top rows are removed
          * and the lower identity matrix is shifted to the left. For negative p the bottom rows are removed and the
          * upper identity is shifted the right
          *
          * @param p Order of the shifted quasi identity
          */
         virtual SparseMatrix shiftId(const int p) const;

         /**
          * @brief Get the derivative of order p
          *
          * @param nBC  Number of boundary
          * @param p    Order of the derivative
          */
         virtual SparseMatrix diff(const int nBC, const int p) const = 0;

         /**
          * @brief Get the quasi inverse derivative for D^-p D^q
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         virtual SparseMatrix qDiff(const int p, const int q) const = 0;
         
      protected:
         /**
          * @brief Get size of spectral space
          */
         int basisN() const;

      private:
         /**
          * @brief Size of the polynomial basis
          */
         int mBasisN;
   };

}
}

#endif // IOPERATOR_HPP
