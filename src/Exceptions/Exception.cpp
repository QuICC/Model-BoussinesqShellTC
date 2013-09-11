/** 
 * @file Exception.cpp
 * @brief Definitions of the Exception methods.
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
