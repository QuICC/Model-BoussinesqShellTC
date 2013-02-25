/** \file TransformSelector.hpp
 *  \brief Definition of some useful typedefs for the variables used in the code
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
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator3D.hpp"

namespace GeoMHDiSCC {

   namespace Transform {

      template<Dimensions::Transform::Id TId> struct TransformSelector;

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
   }

   namespace Parallel {
      typedef Communicator3D<Datatypes::ScalarSelector<Dimensions::Transform::TRA1D>::FwdType, Datatypes::ScalarSelector<Dimensions::Transform::TRA1D>::BwdType, Datatypes::ScalarSelector<Dimensions::Transform::TRA2D>::FwdType, Datatypes::ScalarSelector<Dimensions::Transform::TRA2D>::BwdType, Datatypes::ScalarSelector<Dimensions::Transform::TRA3D>::FwdType, Datatypes::ScalarSelector<Dimensions::Transform::TRA3D>::BwdType> CommunicatorType;
   }

   namespace Transform {

      /// Typedef for a TransformCoordinatorType
      typedef Transform3DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA1D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;

   }
}

#endif // TRANSFORMSELECTOR_HPP
