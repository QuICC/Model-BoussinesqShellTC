/** 
 * @file UnitOperator.hpp
 * @brief Simple unit spectral operator
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

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Simple unit spectral operator
    */
   class UnitOperator
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
         virtual ~UnitOperator();

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

#endif // UNITOPERATOR_HPP
