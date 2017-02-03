/** 
 * @file CflTags.hpp
 * @brief Definitions and names use by the CFL writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CFLTAGS_HPP
#define CFLTAGS_HPP

// System includes
//
#include <string>

// External includes
//

// Project includes
//

namespace QuICC {

namespace IoAscii {

   /**
    * @brief Definitions and names use by the CFL writer
    */
   class CflTags
   {
      public:
         /**
          * @brief HEADER part for Cfl file
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for Cfl file
          */
         static const std::string   VERSION;

         /**
          * @brief NAME of Cfl file
          */
         static const std::string   NAME;

         /**
          * @brief TYPE of Cfl file
          */
         static const std::string   TYPE;

         /**
          * @brief EXTENSION of Cfl file
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         CflTags();

         /**
         * @brief Destructor
         */
         ~CflTags();
   };
}
}

#endif // CFLTAGS_HPP
