/** 
 * @file StatisticsPart.hpp 
 * @brief Implementation of the statistics part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATISTICSPART_HPP
#define STATISTICSPART_HPP

// Configuration includes
//

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "IoConfig/ConfigParts/IConfigurationPart.hpp"

namespace QuICC {

namespace IoConfig {

   /**
    * @brief Implementation of the statistics part of the configuration file
    */
   class StatisticsPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         StatisticsPart();

         /**
          * @brief Destructor
          */
         virtual ~StatisticsPart();

         /**
          * @brief Check compatibility of data
          *
          * \mhdBug Check is not yet implemented
          */
         virtual void checkData();
         
      protected:
         /**
          * @brief Tag name of the parent node
          */
         static const std::string  PARENTTAG;

         /**
          * @brief Initialise component
          */
         void init();

      private:
   };

   /// Typedef for a shared pointer of a IO part
   typedef SharedPtrMacro<StatisticsPart> SharedStatisticsPart;

   /// Typedef for a const shared pointer of a IO part
   typedef SharedPtrMacro<const StatisticsPart> SharedCStatisticsPart;

}
}

#endif // STATISTICSPART_HPP
