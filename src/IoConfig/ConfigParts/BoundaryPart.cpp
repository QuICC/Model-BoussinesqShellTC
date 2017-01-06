/** 
 * @file BoundaryPart.cpp
 * @brief Source of the implementation of the boundary part of the configuration
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
#include "IoConfig/ConfigParts/BoundaryPart.hpp"

// Project includes
//

namespace QuICC {

namespace IoConfig {

   const std::string BoundaryPart::PARENTTAG = "boundary";

   BoundaryPart::BoundaryPart(const std::vector<std::string>& names)
      : IConfigurationPart(BoundaryPart::PARENTTAG)
   {
      this->init(names);
   }

   BoundaryPart::~BoundaryPart()
   {
   }

   void BoundaryPart::init(const std::vector<std::string>& names)
   {
      // Get iterator over vector
      std::vector<std::string>::const_iterator  it;
      for(it = names.begin(); it != names.end(); it++)
      {
         this->addIntegerTag(*it, -1);
      }
   }

   void BoundaryPart::checkData()
   {
   }

}
}
