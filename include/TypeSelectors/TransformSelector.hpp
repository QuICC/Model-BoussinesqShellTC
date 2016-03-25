/** 
 * @file TransformSelector.hpp
 * @brief Typedefs to setup the correct transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMSELECTOR_HPP
#define TRANSFORMSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/FftSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      template<Dimensions::Transform::Id TId> struct TransformSelector;

   }
}

   // Configure code to use TTT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TTT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TTT

   // Configure code to use TFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TFT

   // Configure code to use TFF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TFF

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF

      // Configure code to use FFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_FFF

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_FFF

   // Configure code to use CFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_CFT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::CylinderChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_CFT

   // Configure code to use AFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_AFT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::AnnulusChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_AFT

   // Configure code to use BLF scheme
   #if defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_BLFM

      #if defined GEOMHDISCC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined GEOMHDISCC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined GEOMHDISCC_ALEGTRA_MATRIX

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::SphereChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined GEOMHDISCC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined GEOMHDISCC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined GEOMHDISCC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_BLFM

   // Configure code to use SLF scheme
   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM

      #if defined GEOMHDISCC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined GEOMHDISCC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined GEOMHDISCC_ALEGTRA_MATRIX

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ShellChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined GEOMHDISCC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined GEOMHDISCC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined GEOMHDISCC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM

   // Configure code to use WFT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WFT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef CylinderWorlandTransform Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WFT

   // Configure code to use WLF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_WLF

      #if defined GEOMHDISCC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined GEOMHDISCC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined GEOMHDISCC_ALEGTRA_MATRIX

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef SphereWorlandTransform Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined GEOMHDISCC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined GEOMHDISCC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined GEOMHDISCC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WLF

   // Configure code to use TT scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TT

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::ChebyshevSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TT

   // Configure code to use TF scheme
   #ifdef GEOMHDISCC_SPATIALSCHEME_TF

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_TF

#endif // TRANSFORMSELECTOR_HPP
