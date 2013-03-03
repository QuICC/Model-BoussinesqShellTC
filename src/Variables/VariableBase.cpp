/** \file VariableBase.cpp
 *  \brief Base of the implementation of the variables
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Variables/VariableBase.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Datatypes {

   VariableBase::VariableBase(SharedResolution spRes)
      : mspRes(spRes)
   {
   }

   VariableBase::~VariableBase()
   {
   }

}
}
