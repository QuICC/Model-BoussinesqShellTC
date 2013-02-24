/** \file TimesteppingPart.cpp
 *  \brief Source of the implementation of the parallel part of the configuration
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoConfig/ConfigParts/TimesteppingPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string TimesteppingPart::PARENTTAG = "timestepping";

   TimesteppingPart::TimesteppingPart()
      : ConfigurationPart(TimesteppingPart::PARENTTAG)
   {
      this->init();
   }

   void TimesteppingPart::init()
   {
      this->addFloatTag("time", -1.0);
      this->addFloatTag("timestep", -1.0);
   }

   void TimesteppingPart::checkData()
   {
   }

}
}
