/** \file BoundaryPart.hpp 
 *  \brief Implementation of the boundary part of the configuration file
 *
 *  \mhdBug Needs test
 */

#ifndef BOUNDARYPART_HPP
#define BOUNDARYPART_HPP

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
    * @brief Implementation of the boundary part of the configuration file
    */
   class BoundaryPart: public IConfigurationPart
   {
      public:
         /**
          * @brief Constructor
          *
          * @param names Names of the parameters
          */
         BoundaryPart(const std::vector<std::string>& names);

         /**
          * @brief Destructor
          */
         virtual ~BoundaryPart();

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
         void init(const std::vector<std::string>& names);

      private:
   };

   /// Typedef for a shared pointer of a boundary part
   typedef SharedPtrMacro<BoundaryPart> SharedBoundaryPart;

   /// Typedef for a const shared pointer of a boundary part
   typedef SharedPtrMacro<const BoundaryPart> SharedCBoundaryPart;

}
}

#endif // BOUNDARYPART_HPP
