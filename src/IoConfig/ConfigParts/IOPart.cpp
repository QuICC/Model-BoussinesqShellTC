/** \file IOPart.cpp
 *  \brief Source of the implementation of the IO part of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoConfig/ConfigParts/IOPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string IOPart::PARENTTAG = "io";

   IOPart::IOPart()
      : IConfigurationPart(IOPart::PARENTTAG)
   {
      this->init();
   }

   void IOPart::init()
   {
      this->addFloatTag("arate", 0.0);

      this->addFloatTag("srate", 0.0);
   }

   void IOPart::checkData()
   {
   }

}
}
