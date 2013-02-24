/** \file ParallelPart.cpp
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
#include "IoConfig/ConfigParts/ParallelPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string ParallelPart::PARENTTAG = "parallel";

   ParallelPart::ParallelPart()
      : ConfigurationPart(ParallelPart::PARENTTAG)
   {
      this->init();
   }

   void ParallelPart::init()
   {
      this->addIntegerTag("cpus", -1);
   }

   void ParallelPart::checkData()
   {
   }

}
}
