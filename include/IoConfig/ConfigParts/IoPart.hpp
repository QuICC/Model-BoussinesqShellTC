/** 
 * @file IoPart.hpp 
 * @brief Implementation of the IO part of the configuration file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef IOPART_HPP
#define IOPART_HPP

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
    * @brief Implementation of the IO part of the configuration file
    */
   class IoPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         IoPart();

         /**
          * @brief Destructor
          */
         virtual ~IoPart();

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
   typedef SharedPtrMacro<IoPart> SharedIoPart;

   /// Typedef for a const shared pointer of a IO part
   typedef SharedPtrMacro<const IoPart> SharedCIoPart;

}
}

#endif // IOPART_HPP
