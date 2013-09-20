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
#include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      template<Dimensions::Transform::Id TId> struct TransformSelector;

      // Configure code to use TTT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TTT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef ChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef ChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TTT

      // Configure code to use TFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef ChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFT

      // Configure code to use TFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef ChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF

      // Configure code to use FFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_FFF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_FFF

      // Configure code to use CFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_CFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef CylindricalChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_CFT

      // Configure code to use SLF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_SLF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef SphericalChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef AssociatedLegendreTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_SLF

      // Configure code to use WFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef CylindricalChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_WFT

      // Configure code to use WLF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WLF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef SphericalChebyshevTransformType Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef AssociatedLegendreTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftTransformType Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_WLF
   }

   namespace Parallel {
      typedef Communicator<Dimensions::THREED, Datatypes::ScalarSelector> CommunicatorType;
   }

   namespace Transform {

      /// Typedef for a TransformCoordinatorType
      typedef Transform3DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, TransformSelector<Dimensions::Transform::TRA3D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;

   }
}

#endif // TRANSFORMSELECTOR_HPP
