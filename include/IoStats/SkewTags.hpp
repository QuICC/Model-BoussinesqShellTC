/** 
 * @file SkewTags.hpp
 * @brief Definitions and names use by the average writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SKEWTAGS_HPP
#define SKEWTAGS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace IoStats {

   /**
    * @brief Definitions and names use by the energy writers
    */
   class SkewTags
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
         SkewTags();

         /**
         * @brief Destructor
         */
         ~SkewTags();
   };
}
}

#endif // SKEWTAGS_HPP
