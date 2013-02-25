/** \file BoundaryPart.cpp
 *  \brief Source of the implementation of the boundary part of the configuration
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

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string BoundaryPart::PARENTTAG = "boundary";

   BoundaryPart::BoundaryPart()
      : IConfigurationPart(BoundaryPart::PARENTTAG)
   {
      this->init();
   }

   BoundaryPart::~BoundaryPart()
   {
   }

   void BoundaryPart::init()
   {
      this->addIntegerTag("codensity", 0);

      this->addIntegerTag("velocity", 0);

      this->addIntegerTag("vorticity", 0);

      this->addIntegerTag("magnetic", 0);

      this->addIntegerTag("helicity", 0);
   }

   void BoundaryPart::checkData()
   {
   }

}
}
