/** 
 * @file TransformSelector.hpp
 * @brief Typedefs to setup the correct transforms
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "FastTransforms/FftwTransform.hpp"
#include "FastTransforms/ChebyshevFftwTransform.hpp"
#include "FastTransforms/CylindricalChebyshevFftwTransform.hpp"
#include "FastTransforms/SphericalChebyshevFftwTransform.hpp"
#include "PolynomialTransforms/AssociatedLegendreTransform.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      template<Dimensions::Transform::Id TId> struct TransformSelector;

      // Configure code to use TFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef ChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevFftwTransform Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFT

      // Configure code to use TFF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_TFF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef ChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftwTransform Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_TFF

      // Configure code to use CFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_CFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef CylindricalChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevFftwTransform Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_CFT

      // Configure code to use WFT scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WFT
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef CylindricalChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef FftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef ChebyshevFftwTransform Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_WFT

      // Configure code to use SLF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_SLF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef SphericalChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef AssociatedLegendreTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftwTransform Type;
         };
      #endif //GEOMHDISCC_SPATIALSCHEME_SLF

      // Configure code to use WLF scheme
      #ifdef GEOMHDISCC_SPATIALSCHEME_WLF
         template<> struct TransformSelector<Dimensions::Transform::TRA1D>
         {
            /// Typedef for the first transform
            typedef SphericalChebyshevFftwTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA2D>
         {
            /// Typedef for the second transform
            typedef AssociatedLegendreTransform Type;
         };

         template<> struct TransformSelector<Dimensions::Transform::TRA3D>
         {
            /// Typedef for the third transform
            typedef FftwTransform Type;
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
