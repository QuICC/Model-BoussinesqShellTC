/** 
 * @file RunPart.cpp
 * @brief Source of the implementation of the run part of the configuration
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
#include "IoConfig/ConfigParts/RunPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string RunPart::PARENTTAG = "run";

   RunPart::RunPart()
      : IConfigurationPart(RunPart::PARENTTAG)
   {
      this->init();
   }

   RunPart::~RunPart()
   {
   }

   void RunPart::init()
   {
      this->addFloatTag("sim", 0.0);
      this->addFloatTag("wall", 0.0);
   }

   void RunPart::checkData()
   {
   }

}
}
