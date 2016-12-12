/** 
 * @file TransformCommSelector.hpp
 * @brief Typedefs to setup the correct transforms communicator and coordinator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMCOMMSELECTOR_HPP
#define TRANSFORMCOMMSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformCoordinators/Transform3DCoordinator.hpp"
#else
   #include "TransformCoordinators/Transform2DCoordinator.hpp"
#endif //GEOMHDISCC_SPATIALDIMENSION_3D
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator.hpp"

namespace GeoMHDiSCC {

   namespace Parallel {
      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         typedef Communicator<Dimensions::THREED, Datatypes::ScalarSelector> CommunicatorType;
      #else
         typedef Communicator<Dimensions::TWOD, Datatypes::ScalarSelector> CommunicatorType;
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D
   }

   namespace Transform {
      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         /// Typedef for a TransformCoordinatorType
         typedef Transform3DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, TransformSelector<Dimensions::Transform::TRA3D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;
      #else
         /// Typedef for a TransformCoordinatorType
         typedef Transform2DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D
   }
}

#endif // TRANSFORMCOMMSELECTOR_HPP
