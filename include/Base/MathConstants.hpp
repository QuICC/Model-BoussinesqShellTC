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

namespace QuICC {

namespace Math {

   /**
    * @brief The constant \f$\pi\f$
    */
   const MHDFloat PI = std::acos(-1);

   /**
    * @brief Pure imaginary value I
    */
   const MHDComplex cI = MHDComplex(0.0, 1.0);

}
}

#endif // MATHCONSTANTS_HPP
