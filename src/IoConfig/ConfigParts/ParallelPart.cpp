/** 
 * @file ParallelPart.cpp
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
#include "IoConfig/ConfigParts/ParallelPart.hpp"

// Project includes
//

namespace QuICC {

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
