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
         // Chebyshev ordering r = (1, 0)
         //int k_ = k;
         // Reverse Chebyshev grid r = (0, 1)
         int k_ = size-k-1;

         igrid(k_) = precision::cos((Precision::PI)*(internal::MHDFloat(k)+ MHD_MP(0.5))/internal::MHDFloat(size));

         igrid(k_) = precision::sqrt((igrid(k_) + MHD_MP(1.0))/MHD_MP(2.0));
      }

      iweights.setConstant(Precision::PI/(MHD_MP(2.0)*internal::MHDFloat(size)));

      // Copy internal precision values into input arrays
      grid = Precision::cast(igrid);
      weights = Precision::cast(iweights);
   }

}
