/** 
 * @file IoPart.cpp
 * @brief Source of the implementation of the IO part of the configuration
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
#include "IoConfig/ConfigParts/IoPart.hpp"

// Project includes
//

namespace QuICC {

namespace IoConfig {

   const std::string IoPart::PARENTTAG = "io";

   IoPart::IoPart()
      : IConfigurationPart(IoPart::PARENTTAG)
   {
      this->init();
   }

   IoPart::~IoPart()
   {
   }

   void IoPart::init()
   {
      this->addFloatTag("ascii", 0.0);

      this->addFloatTag("hdf5", 0.0);
   }

   void IoPart::checkData()
   {
   }

}
}
