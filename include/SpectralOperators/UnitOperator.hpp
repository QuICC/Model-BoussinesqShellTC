/** 
 * @file UnitOperator.hpp
 * @brief Implementation of the unit spectral operator with quasi inverses
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef UNITOPERATOR_HPP
#define UNITOPERATOR_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IOperator.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the unit spectral operator with quasi inverses
    */
   class UnitOperator: IOperator
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         UnitOperator(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~UnitOperator();

         /**
          * @brief Get the identity matrix
          *
          * @param p Order of the quasi identity
          */
         SparseMatrix id(const int p) const;

         /**
          * @brief Get the derivative of order p
          *
          * @param nBC  Number of boundary
          * @param p Order of the derivative
          */
         SparseMatrix diff(const int nBC, const int p) const;

         /**
          * @brief Get the quasi inverse derivative of order p
          *
          * @param p    Order of the quasi inverse
          * @param q    Order of the derivative
          */
         SparseMatrix qDiff(const int p, const int q) const;
         
      protected:

      private:
   };

}
}

#endif // UNITOPERATOR_HPP
