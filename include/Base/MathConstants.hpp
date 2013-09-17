/**
 * @file MathConstants.hpp
 * @brief Definition of some useful math constants
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MATHCONSTANTS_HPP
#define MATHCONSTANTS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Contains some useful math constants
    */
   class MathConstants
   {
      public:
         /**
          * @brief The constant \f$\pi\f$
          */
         static const MHDFloat PI;

         /**
          * @brief Pure imaginary value I
          */
         static const MHDComplex cI;
         
      protected:

      private:
         /**
          * @brief Empty constructor
          */
         MathConstants();

         /**
          * @brief Empty Destructor
          */
         ~MathConstants();
   };

}

#endif // MATHCONSTANTS_HPP
