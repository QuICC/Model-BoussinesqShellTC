/** 
 * @file AngularMomentumTags.hpp
 * @brief Definitions and names use by the angular momentum writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ANGULARMOMENTUMTAGS_HPP
#define ANGULARMOMENTUMTAGS_HPP

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
    * @brief Definitions and names use by the angular momentum writers
    */
   class AngularMomentumTags
   {
      public:
         /**
          * @brief HEADER part
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         AngularMomentumTags();

         /**
         * @brief Destructor
         */
         ~AngularMomentumTags();
   };
}
}

#endif // ANGULARMOMENTUMTAGS_HPP
