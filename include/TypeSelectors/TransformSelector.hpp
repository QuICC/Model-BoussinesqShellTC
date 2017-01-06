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

namespace QuICC {

   namespace Transform {

      template<Dimensions::Transform::Id TId> struct TransformSelector;

   }
}

   // Configure code to use TTT scheme
   #ifdef QUICC_SPATIALSCHEME_TTT

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_TTT

   // Configure code to use TFT scheme
   #ifdef QUICC_SPATIALSCHEME_TFT

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_TFT

   // Configure code to use TFF scheme
   #ifdef QUICC_SPATIALSCHEME_TFF

      namespace QuICC {
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
      #endif //QUICC_SPATIALSCHEME_TFF

      // Configure code to use FFF scheme
      #ifdef QUICC_SPATIALSCHEME_FFF

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_FFF

   // Configure code to use CFT scheme
   #ifdef QUICC_SPATIALSCHEME_CFT

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_CFT

   // Configure code to use AFT scheme
   #ifdef QUICC_SPATIALSCHEME_AFT

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_AFT

   // Configure code to use BLF scheme
   #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM

      #if defined QUICC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined QUICC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined QUICC_ALEGTRA_MATRIX

      namespace QuICC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::SphereChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined QUICC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined QUICC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined QUICC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_BLFM

   // Configure code to use SLF scheme
   #if defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM

      #if defined QUICC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined QUICC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined QUICC_ALEGTRA_MATRIX

      namespace QuICC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::ShellChebyshevSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined QUICC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined QUICC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined QUICC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //defined QUICC_SPATIALSCHEME_SLFL || defined QUICC_SPATIALSCHEME_SLFM

   // Configure code to use WFT scheme
   #ifdef QUICC_SPATIALSCHEME_WFT

      #if defined QUICC_WORLANDTRA_MATRIX
         #include "PolynomialTransforms/CylinderWorlandTransform.hpp"
      #elif defined QUICC_WORLANDTRA_FLY
         #include "PolynomialTransforms/CylinderWorlandFlyTransform.hpp"
      #endif //defined QUICC_WORLANDTRA_MATRIX

      namespace QuICC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               #if defined QUICC_WORLANDTRA_MATRIX
                  typedef CylinderWorlandTransform Type;
               #elif defined QUICC_WORLANDTRA_FLY 
                  typedef CylinderWorlandFlyTransform Type;
               #endif //defined QUICC_WORLANDTRA_MATRIX
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
   #endif //QUICC_SPATIALSCHEME_WFT

   // Configure code to use WLF scheme
   #if defined QUICC_SPATIALSCHEME_WLFL || defined QUICC_SPATIALSCHEME_WLFM

      #if defined QUICC_WORLANDTRA_MATRIX
         #include "PolynomialTransforms/SphereWorlandTransform.hpp"
      #elif defined QUICC_WORLANDTRA_FLY
         #include "PolynomialTransforms/SphereWorlandFlyTransform.hpp"
      #endif //defined QUICC_WORLANDTRA_MATRIX

      #if defined QUICC_ALEGTRA_MATRIX
         #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
      #elif defined QUICC_ALEGTRA_FLY 
         #include "PolynomialTransforms/AssociatedLegendreFlyTransform.hpp"
      #endif //defined QUICC_ALEGTRA_MATRIX

      namespace QuICC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               #if defined QUICC_WORLANDTRA_MATRIX
                  typedef SphereWorlandTransform Type;
               #elif defined QUICC_WORLANDTRA_FLY 
                  typedef SphereWorlandFlyTransform Type;
               #endif //defined QUICC_WORLANDTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               #if defined QUICC_ALEGTRA_MATRIX
                  typedef AssociatedLegendreTransform Type;
               #elif defined QUICC_ALEGTRA_FLY 
                  typedef AssociatedLegendreFlyTransform Type;
               #endif //defined QUICC_ALEGTRA_MATRIX
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //defined QUICC_SPATIALSCHEME_WLFL || defined QUICC_SPATIALSCHEME_WLFM

   // Configure code to use TT scheme
   #ifdef QUICC_SPATIALSCHEME_TT

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_TT

   // Configure code to use TF scheme
   #ifdef QUICC_SPATIALSCHEME_TF

      namespace QuICC {
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
   #endif //QUICC_SPATIALSCHEME_TF

#endif // TRANSFORMSELECTOR_HPP
