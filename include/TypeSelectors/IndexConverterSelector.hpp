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

namespace GeoMHDiSCC {

   namespace Parallel {

      /**
       * @brief Template class to selec the index converter for different transform
       */
      template<Dimensions::Transform::Id TId> struct IndexConverterSelector;

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

      /**
       * @brief Specialialised IndexConverterSelector for the second transform
       */
      // Configure index converter on spherical harmonics basis with l spectral ordering
   #if defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_WLFL
      #include "Communicators/Converters/SHlIndexConv.hpp"

      namespace GeoMHDiSCC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef SHlIndexConv  Type;
            };
         }
      }
      // Configure index converter on spherical harmonics basis with m spectral ordering
   #elif defined GEOMHDISCC_SPATIALSCHEME_BLFM  || defined GEOMHDISCC_SPATIALSCHEME_SLFM || defined GEOMHDISCC_SPATIALSCHEME_WLFM
      #include "Communicators/Converters/SHmIndexConv.hpp"

      namespace GeoMHDiSCC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef SHmIndexConv  Type;
            };
         }
      }
      // Configure index converter on in-order C2C FFT values (plus-minus frequency order)
   #elif defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_FFF

      #include "Communicators/Converters/PMIndexConv.hpp"

      namespace GeoMHDiSCC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef PMIndexConv  Type;
            };
         }
      }
   #else

      namespace GeoMHDiSCC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef NoIndexConv  Type;
            };
         }
      }
   #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_WLFL

#endif // INDEXCONVERTERSELCTOR_HPP
