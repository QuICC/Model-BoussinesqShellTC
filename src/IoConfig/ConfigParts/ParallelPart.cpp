/** 
 * @file ParallelPart.cpp
 * @brief Source of the implementation of the parallel part of the configuration
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
#include "IoConfig/ConfigParts/ParallelPart.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string ParallelPart::PARENTTAG = "parallel";

   ParallelPart::ParallelPart()
      : IConfigurationPart(ParallelPart::PARENTTAG)
   {
      this->init();
   }

   ParallelPart::~ParallelPart()
   {
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
