/** 
 * @file TimesteppingPart.cpp
 * @brief Source of the implementation of the parallel part of the configuration
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
#include "IoConfig/ConfigParts/TimesteppingPart.hpp"

// Project includes
//

namespace QuICC {

namespace IoConfig {

   const std::string TimesteppingPart::PARENTTAG = "timestepping";

   TimesteppingPart::TimesteppingPart()
      : IConfigurationPart(TimesteppingPart::PARENTTAG)
   {
      this->init();
   }

   TimesteppingPart::~TimesteppingPart()
   {
   }

   void TimesteppingPart::init()
   {
      this->addFloatTag("time", -1.0);
      this->addFloatTag("timestep", -1.0);
      this->addFloatTag("error", -1.0);
   }

   void TimesteppingPart::checkData()
   {
   }

}
}
