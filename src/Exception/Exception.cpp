/** \file Exception.cpp
 *  \brief Definitions of the Exception methods.
 */

// System includes
//

// External includes
//

// Class include
//
#include "Exception/Exception.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   Exception::Exception(const std::string& msg)
      : std::runtime_error(msg)   
   {
   }

}
