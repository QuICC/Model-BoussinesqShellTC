/** 
 * @file EnstrophySpectrumTags.hpp
 * @brief Definitions and names use by the enstrophy spectrum writers
 */

#ifndef QUICC_IO_VARIABLE_ENSTROPHYSPECTRUMTAGS_HPP
#define QUICC_IO_VARIABLE_ENSTROPHYSPECTRUMTAGS_HPP

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
    * @brief Definitions and names use by the enstrophy spectrum writers
    */
   class EnstrophySpectrumTags
   {
      public:
         /**
          * @brief HEADER part for enstrophy spectrum files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for enstrophy spectrum files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of enstrophy spectrum files
          */
         static const std::string   BASENAME;

         /**
          * @brief BASENAME of L enstrophy spectrum files
          */
         static const std::string   LBASENAME;

         /**
          * @brief BASENAME of M enstrophy spectrum files
          */
         static const std::string   MBASENAME;

         /**
          * @brief EXTENSION of enstrophy spectrum files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         EnstrophySpectrumTags();

         /**
         * @brief Destructor
         */
         ~EnstrophySpectrumTags();
   };
}
}

#endif // QUICC_IO_VARIABLE_ENSTROPHYSPECTRUMTAGS_HPP
