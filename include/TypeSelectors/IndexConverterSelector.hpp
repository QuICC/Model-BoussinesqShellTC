/** 
 * @file IndexConverterSelector.hpp
 * @brief Preprocessor macros to setup the index converters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef INDEXCONVERTERSELCTOR_HPP
#define INDEXCONVERTERSELCTOR_HPP

//
// Create index converter macro macros
//

// include serial communicators
#include "Communicators/Converters/NoIndexConv.hpp"
#include "Communicators/Converters/SHIndexConv.hpp"

namespace GeoMHDiSCC {

   namespace Parallel {

      /**
       * @brief Template class to selec the index converter for different transform
       */
      template<Dimensions::Transform::Id TId> struct IndexConverterSelector;

      /**
       * @brief Specialialised IndexConverterSelector for the second transform
       */
      // Configure code to use SLF scheme
      #if defined GEOMHDISCC_SPATIALSCHEME_SLF || defined GEOMHDISCC_SPATIALSCHEME_WLF
         template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the type of the index converter
            typedef SHIndexConv  Type;
         };
      #else //GEOMHDISCC_SPATIALSCHEME_SLF || GEOMHDISCC_SPATIALSCHEME_WLF
         template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the type of the index converter
            typedef NoIndexConv  Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_SLF || GEOMHDISCC_SPATIALSCHEME_WLF


      /**
       * @brief Specialialised IndexConverterSelector for the second transform
       */
      template <> struct IndexConverterSelector<Dimensions::Transform::TRA3D>
      {
         /// Typedef for the type of the index converter
         typedef NoIndexConv  Type;
      };
   }
}

#endif // INDEXCONVERTERSELCTOR_HPP
