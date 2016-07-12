/** 
 * @file WorlandChebyshevRule.hpp
 * @brief Implementation of the Worland Chebyshev quadrature rule
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef WORLANDCHEBYSHEVRULE_HPP
#define WORLANDCHEBYSHEVRULE_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Precision.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the Worland Chebyshev quadrature rule
    */
   class WorlandChebyshevRule
   {
      public:
         /**
          * @brief Compute the quadrature and return standard precision values
          */
         static void computeQuadrature(Array& grid, Array& weights, const int size);

         /**
          * @brief Compute the quadrature and return standard and internal precision values
          */
         static void computeQuadrature(Array& grid, Array& weights, internal::Array& igrid, internal::Array& iweights, const int size);
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         WorlandChebyshevRule();

         /**
          * @brief Empty destructor
          */
         virtual ~WorlandChebyshevRule();
   };

}

#endif // WORLANDCHEBYSHEVRULE_HPP
