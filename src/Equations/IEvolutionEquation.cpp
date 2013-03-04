/** \file IEvolutionEquation.cpp
 *  \brief Source of building block for the implementation of a time dependend evolution equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/IEvolutionEquation.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Equations {

   IEvolutionEquation::IEvolutionEquation(SharedIEquationParameters spEqParams)
      : EquationData(spEqParams)
   {
   }

   IEvolutionEquation::~IEvolutionEquation()
   {
   }

   void IEvolutionEquation::init()
   {
      this->setCoupling();
   }
}
}
