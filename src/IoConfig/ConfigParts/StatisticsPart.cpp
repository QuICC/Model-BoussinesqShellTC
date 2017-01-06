/** 
 * @file StatisticsPart.cpp
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
#include "IoConfig/ConfigParts/StatisticsPart.hpp"

// Project includes
//

namespace QuICC {

namespace IoConfig {

   const std::string StatisticsPart::PARENTTAG = "statistics";

   StatisticsPart::StatisticsPart()
      : IConfigurationPart(StatisticsPart::PARENTTAG)
   {
      this->init();
   }

   StatisticsPart::~StatisticsPart()
   {
   }

   void StatisticsPart::init()
   {
      this->addFloatTag("output_rate", 0.0);

      this->addFloatTag("time_avg_rate", 0.0);
   }

   void StatisticsPart::checkData()
   {
   }

}
}
