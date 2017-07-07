/** 
 * @file EnergyTags.hpp
 * @brief Definitions and names use by the energy writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef OBSERVABLETAGS_HPP
#define OBSERVABLETAGS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Definitions and names use by the energy writers
    */
   class ObservableTags
   {
      public:
         /**
          * @brief HEADER part for energy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of energy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         ObservableTags();

         /**
         * @brief Destructor
         */
         ~ObservableTags();
   };
}
}

#endif // OBSERVABLETAGS_HPP
