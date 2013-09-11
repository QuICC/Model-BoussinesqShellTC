/** 
 * @file TimesteppingPart.hpp 
 * @brief Implementation of the timestepping part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef TIMESTEPPINGPART_HPP
#define TIMESTEPPINGPART_HPP

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

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Implementation of the timestepping part of the configuration file
    */
   class TimesteppingPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         TimesteppingPart();

         /**
          * @brief Destructor
          */
         virtual ~TimesteppingPart();

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

   /// Typedef for a shared pointer of a timestepping part
   typedef SharedPtrMacro<TimesteppingPart> SharedTimesteppingPart;

   /// Typedef for a const shared pointer of a timestepping part
   typedef SharedPtrMacro<const TimesteppingPart> SharedCTimesteppingPart;

}
}

#endif // TIMESTEPPINGPART_HPP
