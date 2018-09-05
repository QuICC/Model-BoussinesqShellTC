/** 
 * @file Version.cpp
 * @brief Implementation of version class
 * @author Philippe Marti \<philippe.marti@env.ethz.ch\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Base/Version.hpp"

// Project includes
//

namespace QuICC {

   const int Version::MAJOR = QUICC_VERSION_MAJOR;

   const int Version::MINOR = QUICC_VERSION_MINOR;

   const int Version::PATCH = QUICC_VERSION_PATCH;

   std::string Version::version()
   {
      std::stringstream oss;

      oss << Version::MAJOR << "." << Version::MINOR << "." << Version::PATCH;

      return oss.str();
   }

}
