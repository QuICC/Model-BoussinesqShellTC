/** 
 * @file SpectrumTags.hpp
 * @brief Definitions and names use by the energy spectrum writers
 */

#ifndef QUICC_IO_VARIABLE_SPECTRUMTAGS_HPP
#define QUICC_IO_VARIABLE_SPECTRUMTAGS_HPP

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
    * @brief Definitions and names use by the energy spectrum writers
    */
   class SpectrumTags
   {
      public:
         /**
          * @brief HEADER part for energy spectrum files
          */
         static const std::string   HEADER;

         /**
          * @brief VERSION part for energy spectrum files
          */
         static const std::string   VERSION;

         /**
          * @brief BASENAME of energy spectrum files
          */
         static const std::string   BASENAME;

         /**
          * @brief BASENAME of L energy spectrum files
          */
         static const std::string   LBASENAME;

         /**
          * @brief BASENAME of M energy spectrum files
          */
         static const std::string   MBASENAME;

         /**
          * @brief EXTENSION of energy spectrum files
          */
         static const std::string   EXTENSION;

      private:
         /**
         * @brief Empty destructor
         */
         SpectrumTags();

         /**
         * @brief Destructor
         */
         ~SpectrumTags();
   };
}
}

#endif // QUICC_IO_VARIABLE_SPECTRUMTAGS_HPP
