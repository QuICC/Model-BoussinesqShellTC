/**
 * @file Version.hpp
 * @brief Definition of version information
 * @author Philippe Marti \<philippe.marti@env.ethz.ch\>
 */

#ifndef VERSION_HPP
#define VERSION_HPP

#define QUICC_VERSION_MAJOR 0
#define QUICC_VERSION_MINOR 6
#define QUICC_VERSION_PATCH 0

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace QuICC {

   /**
    * @brief Small class storing version information
    */
   class Version 
   {
      public:
         /**
          * @brief Major version number, ie MAJOR.MINOR.PATCH
          */
         static const int MAJOR;

         /**
          * @brief Minor version number, ie MAJOR.MINOR.PATCH
          */
         static const int MINOR;

         /**
          * @brief Patch version number, ie MAJOR.MINOR.PATCH
          */
         static const int PATCH;

         /**
          * @brief Get version number string, ie MAJOR.MINOR.PATCH
          */
         static std::string version();
   };

}

#endif // VERSION_HPP
