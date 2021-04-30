/** 
 * @file EnstrophyTags.hpp
 * @brief Definitions and names use by the enstrophy writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ENERGYTAGS_HPP
#define ENERGYTAGS_HPP

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
    * @brief Definitions and names use by the enstrophy writers
    */
   class EnstrophyTags
   {
      public:
         /**
          * @brief HEADER part for enstrophy files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for enstrophy files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of enstrophy files
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of enstrophy files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         EnstrophyTags();

         /**
         * @brief Destructor
         */
         ~EnstrophyTags();
   };
}
}

#endif // ENERGYTAGS_HPP
