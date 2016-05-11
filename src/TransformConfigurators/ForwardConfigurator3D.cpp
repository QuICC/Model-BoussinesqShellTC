/** 
 * @file ForwardConfigurator3D.cpp
 * @brief Source of the implementation of the base forward configurator in 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Debug/DebuggerMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/ForwardConfigurator3D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardConfigurator3D::integrate2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 2D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD2D);

      TransformCoordinatorType::CommunicatorType::Fwd2DType* pInVar;
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverFwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveForward<Dimensions::Transform::TRA2D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideBwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD2DTRA);

      // Compute integration transform of second dimension
      coord.transform2D().integrate(rOutVar.rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA2D>::Type::IntegratorType::Id>(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD2DTRA);

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdFwd(*pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeFwd(*pInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferBackward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD2D);
   }

}
}
