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
#include "SpectralOperators/UnitOperator.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Interface for a general spectral operator implementation (direct and quasi inverses)
    */
   class IOperator: public UnitOperator
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

      private:
   };

}
}

#endif // IOPERATOR_HPP
