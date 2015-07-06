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
#include "TransformCoordinators/Transform3DCoordinator.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Communicators/Communicator.hpp"

namespace GeoMHDiSCC {

   namespace Parallel {
      #if defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT
         typedef Communicator<Dimensions::TWOD, Datatypes::ScalarSelector> CommunicatorType;
      #else
         typedef Communicator<Dimensions::THREED, Datatypes::ScalarSelector> CommunicatorType;
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT
   }

   namespace Transform {
      #if defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT
         /// Typedef for a TransformCoordinatorType
         typedef Transform2DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;
      #else
         /// Typedef for a TransformCoordinatorType
         typedef Transform3DCoordinator<TransformSelector<Dimensions::Transform::TRA1D>::Type, TransformSelector<Dimensions::Transform::TRA2D>::Type, TransformSelector<Dimensions::Transform::TRA3D>::Type, Parallel::CommunicatorType> TransformCoordinatorType;
      #endif //defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_TT
   }
}

#endif // TRANSFORMCOMMSELECTOR_HPP
