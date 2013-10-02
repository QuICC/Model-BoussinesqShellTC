/** 
 * @file IBoundary.hpp
 * @brief Interface for a general spectral boundary operator implementation
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBOUNDARY_HPP
#define IBOUNDARY_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "BoundaryCondition/BoundaryCondition.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Interface for a general spectral boundary operator implementation
    */
   class IBoundary
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         IBoundary(const int basisN);

         /**
          * @brief Empty Destructor
          */
         virtual ~IBoundary();

         /**
          * @brief Get size of spectral space
          */
         int basisN() const;

         /**
          * @brief Reset spectral operators and set new dimension
          *
          * @param basisN   New size of the spectral basis
          */
         void reset(const int basisN);

         /**
          * @brief Get value at boundary
          *
          * @param pt   boundary point
          */
         virtual Array value(Boundary::BCPosition pt) const = 0;

         /**
          * @brief Get first derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array firstDerivative(Boundary::BCPosition pt) const = 0;

         /**
          * @brief Get second derivative at boundary
          *
          * @param pt   boundary point
          */
         virtual Array secondDerivative(Boundary::BCPosition pt) const = 0;
         
      protected:

      private:
         /**
          * @brief Size of the polynomial basis
          */
         int mBasisN;
   };

}
}

#endif // IBOUNDARY_HPP
