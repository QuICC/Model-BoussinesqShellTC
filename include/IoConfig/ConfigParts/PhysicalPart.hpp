/** \file PhysicalPart.hpp 
 *  \brief Implementation of the physical part of the configuration file
 *
 *  \mhdBug Needs test
 */

#ifndef PHYSICALPART_HPP
#define PHYSICALPART_HPP

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
    * @brief Implementation of the physical part of the configuration file
    */
   class PhysicalPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          *
          * @param names Names of the parameters
          */
         explicit PhysicalPart(const std::vector<std::string>& names);

         /**
          * @brief Destructor
          */
         virtual ~PhysicalPart();

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
          *
          * @param names Names of the parameters
          */
         void init(const std::vector<std::string>& names);

      private:
   };

   /// Typedef for a shared pointer of a physical part
   typedef SharedPtrMacro<PhysicalPart> SharedPhysicalPart;

   /// Typedef for a const shared pointer of a physical part
   typedef SharedPtrMacro<const PhysicalPart> SharedCPhysicalPart;

}
}

#endif // PHYSICALPART_HPP
