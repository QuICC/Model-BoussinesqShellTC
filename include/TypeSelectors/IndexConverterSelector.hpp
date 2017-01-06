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

namespace QuICC {

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
   #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_WLFL
      #include "Communicators/Converters/SHlIndexConv.hpp"

      namespace QuICC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef SHlIndexConv  Type;
            };
         }
      }
      // Configure index converter on spherical harmonics basis with m spectral ordering
   #elif defined QUICC_SPATIALSCHEME_BLFM  || defined QUICC_SPATIALSCHEME_SLFM || defined QUICC_SPATIALSCHEME_WLFM
      #include "Communicators/Converters/SHmIndexConv.hpp"

      namespace QuICC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef SHmIndexConv  Type;
            };
         }
      }
      // Configure index converter on in-order C2C FFT values (plus-minus frequency order)
   #elif defined QUICC_SPATIALSCHEME_TFF || defined QUICC_SPATIALSCHEME_FFF

      #include "Communicators/Converters/PMIndexConv.hpp"

      namespace QuICC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef PMIndexConv  Type;
            };
         }
      }
   #else

      namespace QuICC {
         namespace Parallel {

            template <> struct IndexConverterSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the type of the index converter
               typedef NoIndexConv  Type;
            };
         }
      }
   #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_WLFL

#endif // INDEXCONVERTERSELCTOR_HPP
