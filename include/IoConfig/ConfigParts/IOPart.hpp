/** \file IOPart.hpp 
 *  \brief Implementation of the IO part of the configuration file
 *
 *  \mhdBug Needs test
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
   class IOPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          */
         IOPart();

         /**
          * @brief Destructor
          */
         virtual ~IOPart();

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
   typedef SharedPtrMacro<IOPart> SharedIOPart;

   /// Typedef for a const shared pointer of a IO part
   typedef SharedPtrMacro<const IOPart> SharedCIOPart;

}
}

#endif // IOPART_HPP
