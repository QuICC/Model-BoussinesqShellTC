/** 
 * @file DissipationTags.hpp
 * @brief Definitions and names use by the energy writers
 * @author Nicolò Lardelli \<nicolo.lardelli@erdw.ethz.ch\>
 */

#ifndef DISSIPATIONTAGS_HPP
#define DISSIPATIONTAGS_HPP

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
   class DissipationTags
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
         DissipationTags();

         /**
         * @brief Destructor
         */
         ~DissipationTags();
   };
}
}

#endif // DissipationTAGS_HPP
