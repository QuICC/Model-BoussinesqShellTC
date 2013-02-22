/** \file Exception.cpp
 *  \brief Definitions of the Exception methods.
 */

// System includes
//

// External includes
//

// Class include
//
#include "Exceptions/Exception.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   Exception::Exception(const std::string& msg)
      : std::runtime_error(msg)   
   {
   }

   Exception::~Exception() throw()
   {
   }

}
