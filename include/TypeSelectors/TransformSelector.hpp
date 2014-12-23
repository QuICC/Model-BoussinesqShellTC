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
#include "TransformCoordinators/Transform3DCoordinator.hpp"
#include "TypeSelectors/FftSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator.hpp"

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
   #ifdef GEOMHDISCC_SPATIALSCHEME_BLF

      #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"

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
               typedef AssociatedLegendreTransform Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_BLF

   // Configure code to use SLF scheme
   #if defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM

      #include "PolynomialTransforms/AssociatedLegendreTransform.hpp"

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
               typedef AssociatedLegendreTransform Type;
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
               typedef Fft::CylinderWorlandSelector Type;
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

      namespace GeoMHDiSCC {
         namespace Transform {

            template<> struct TransformSelector<Dimensions::Transform::TRA1D>
            {
               /// Typedef for the first transform
               typedef Fft::SpherWorlandSelector Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA2D>
            {
               /// Typedef for the second transform
               typedef AssociatedLegendreTransform Type;
            };

            template<> struct TransformSelector<Dimensions::Transform::TRA3D>
            {
               /// Typedef for the third transform
               typedef Fft::FftSelector Type;
            };
         }
      }
   #endif //GEOMHDISCC_SPATIALSCHEME_WLF

namespace GeoMHDiSCC {

   namespace Parallel {
      typedef Communicator<Dimensions::THREED, Datatypes::ScalarSelector> CommunicatorType;
   }

   namespace Transform {

      /// Typedef for a TransformCoordinatorType
      typedef Transform3DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, TransformSelector<Dimensions::Transform::TRA3D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;

   }
}

#endif // TRANSFORMSELECTOR_HPP
