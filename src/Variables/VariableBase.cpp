/** 
 * @file VariableBase.cpp
 * @brief Base of the implementation of the variables
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
