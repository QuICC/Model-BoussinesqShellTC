/** 
 * @file KineticEnergyTags.hpp
 * @brief Definitions and names use by the kinetic energy writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef KINETICENERGYTAGS_HPP
#define KINETICENERGYTAGS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Definitions and names use by the Nusselt writer
    */
   class KineticEnergyTags
   {
      public:
         /**
          * @brief HEADER part for Nusselt file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Nusselt file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Nusselt file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Nusselt file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         KineticEnergyTags();

         /**
         * @brief Destructor
         */
         ~KineticEnergyTags();
   };
}
}

#endif // KINETICENERGYTAGS_HPP
