/** 
 * @file SphericalChebyshevBoundary.hpp
 * @brief Implementation of the spectral boundary operator for the chebyshev basis for a spherical radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef SPHERICALCHEBYSHEVBOUNDARY_HPP
#define SPHERICALCHEBYSHEVBOUNDARY_HPP

// System includes
//
#include <set>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "SpectralOperators/IBoundary.hpp"

namespace GeoMHDiSCC {

namespace Spectral {

   /**
    * @brief Implementation of the spectral boundary operators for the chebyshev basis for a spherical radius
    */
   class SphericalChebyshevBoundary: public IBoundary
   {
      public:
         /**
          * @brief Constructor
          *
          * @param basisN   Size of the spectral basis
          */
         SphericalChebyshevBoundary(const int basisN);

         /**
          * @brief Empty Destructor
          */
         ~SphericalChebyshevBoundary();

         /**
          * @brief Get value at boundary
          *
          * @param pt   boundary point
          */
         Array value(Position pt) const;

         /**
          * @brief Get first derivative at boundary
          *
          * @param pt   boundary point
          */
         Array firstDerivative(Position pt) const;

         /**
          * @brief Get second derivative at boundary
          *
          * @param pt   boundary point
          */
         Array secondDerivative(Position pt) const;
         
      protected:

      private:
         /**
          * @brief C factor
          *
          * @param n Order of polynial
          */
         MHDFloat c(const int n) const;

         /**
          * @brief Get array of unit values for left boundary
          */
         Array leftUnit() const;

         /**
          * @brief Get array of unit values for right boundary
          */
         Array rightUnit() const;
   };

}
}

#endif // SPHERICALCHEBYSHEVBOUNDARY_HPP
