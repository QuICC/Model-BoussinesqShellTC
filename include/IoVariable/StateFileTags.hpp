/** 
 * @file StateFileTags.hpp
 * @brief Definitions and names use by the state file readers/writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STATEFILETAGS_HPP
#define STATEFILETAGS_HPP

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
    * @brief Definitions and names use by the state file readers/writers
    */
   class StateFileTags
   {
      public:
         /**
          * @brief HEADER part for State file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for State file
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of State file
          */
         static const std::string   BASENAME;

         /**
          * @brief EXTENSION of State file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         StateFileTags();

         /**
         * @brief Destructor
         */
         ~StateFileTags();
   };
}
}

#endif // STATEFILETAGS_HPP
