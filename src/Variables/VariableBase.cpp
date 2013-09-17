/** 
 * @file VariableBase.cpp
 * @brief Base of the implementation of the variables
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
