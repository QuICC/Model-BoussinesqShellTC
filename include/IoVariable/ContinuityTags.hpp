/** 
 * @file ContinuityTags.hpp
 * @brief Definitions and names use by the continuity writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CONTINUITYTAGS_HPP
#define CONTINUITYTAGS_HPP

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
    * @brief Definitions and names use by the continuity writer
    */
   class ContinuityTags
   {
      public:
         /**
          * @brief HEADER part for Continuity file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Continuity file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of Continuity file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of Continuity file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         ContinuityTags();

         /**
         * @brief Destructor
         */
         ~ContinuityTags();
   };
}
}

#endif // CONTINUITYTAGS_HPP
