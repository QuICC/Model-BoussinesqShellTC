/** \file IEquationParameters.cpp
 *  \brief Source of the implementation of the non dimensional parameters
 */

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IEquationParameters.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IEquationParameters::IEquationParameters()
   {
   }

   IEquationParameters::~IEquationParameters()
   {
   }

   MHDFloat IEquationParameters::nd(NonDimensional::Id name) const
   {
      return this->mND.at(name);
   }

}
}
