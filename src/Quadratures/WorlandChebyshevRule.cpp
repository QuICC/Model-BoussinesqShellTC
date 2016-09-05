/** 
 * @file WorlandChebyshevRule.cpp
 * @brief Source of the Worland Chebyshev quadrature
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Quadratures/WorlandChebyshevRule.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   void WorlandChebyshevRule::computeQuadrature(Array& grid, Array& weights, const int size)
   {
      // Internal grid and weights arrays
      internal::Array   igrid;
      internal::Array   iweights;

      WorlandChebyshevRule::computeQuadrature(grid, weights, igrid, iweights, size);
   }

   void WorlandChebyshevRule::computeQuadrature(Array& grid, Array& weights, internal::Array& igrid, internal::Array& iweights, const int size)
   {
      // Internal grid and weights arrays
      igrid.resize(size);
      iweights.resize(size);

      for(int k = 0; k < size; k++)
      {
         igrid(k) = precision::cos((Precision::PI)*(internal::MHDFloat(k)+ MHD_MP(0.5))/internal::MHDFloat(size));

         igrid(k) = precision::sqrt((igrid(k) + MHD_MP(1.0))/MHD_MP(2.0));
      }

      iweights.setConstant(Precision::PI/(MHD_MP(2.0)*MHD_MP(size)));

      // Copy internal precision values into input arrays
      grid = Precision::cast(igrid);
      weights = Precision::cast(iweights);
   }

}
