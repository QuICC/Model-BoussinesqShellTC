/** 
 * @file PhysicalModel.cpp
 * @brief Source of the Boussinesq thermal convection in a spherical shell (Toroidal/Poloidal formulation) without coupled solve (standard implementation)
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "QuICC/Model/Boussinesq/Shell/TC/Explicit/PhysicalModel.hpp"

// Project includes
//

namespace QuICC {

namespace Model {

namespace Boussinesq {

namespace Shell {

namespace TC {

namespace Explicit {

   std::string PhysicalModel::PYMODULE()
   {
      return "boussinesq.shell.tc.explicit.physical_model";
   }

}
}
}
}
}
}
