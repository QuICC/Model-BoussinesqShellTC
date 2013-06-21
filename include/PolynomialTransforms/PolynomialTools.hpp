/** \file PolynomialTools.hpp
 *  \brief Definition of some useful constants and tools for polynomial transforms
 */

#ifndef POLYNOMIALTOOLS_HPP
#define POLYNOMIALTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief Contains some useful constants and tools for polynomial transforms
    */
   class PolynomialTools
   {
      public:
         /**
          * @brief Compute the dealiased size using the 3/2a-rule
          *
          * @param size Size to dealias
          */
         static int dealias(const int size);
         
      protected:

      private:
         /**
          * @brief Standard dealiasing factor (usually 3/2)
          */
         static const MHDFloat STD_DEALIASING;

         /**
          * @brief Empty constructor
          */
         PolynomialTools();

         /**
          * @brief Empty Destructor
          */
         ~PolynomialTools();

   };

}
}

#endif // POLYNOMIALTOOLS_HPP
