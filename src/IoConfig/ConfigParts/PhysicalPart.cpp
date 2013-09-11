/** 
 * @file PhysicalPart.cpp
 * @brief Source of the implementation of the physical part of the configuration
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
#include "IoConfig/ConfigParts/PhysicalPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {
   
namespace IoConfig {

   const std::string PhysicalPart::PARENTTAG = "physical";

   PhysicalPart::PhysicalPart(const std::vector<std::string>& names)
      : IConfigurationPart(PhysicalPart::PARENTTAG)
   {
      this->init(names);
   }

   PhysicalPart::~PhysicalPart()
   {
   }

   void PhysicalPart::init(const std::vector<std::string>& names)
   {
      // Get iterator over vector
      std::vector<std::string>::const_iterator  it;
      for(it = names.begin(); it != names.end(); it++)
      {
         this->addFloatTag(*it, -1.0);
      }
   }

   void PhysicalPart::checkData()
   {
   }

}
}
