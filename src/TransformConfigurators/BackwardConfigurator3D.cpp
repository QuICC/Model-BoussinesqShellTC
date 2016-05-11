/** 
 * @file BackwardConfigurator3D.cpp
 * @brief Source of the implementation of the base backward configurator in 3D space
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
#include "TransformConfigurators/BackwardConfigurator3D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void BackwardConfigurator3D::project2D(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Project 2D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd2DType* pInVar;
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveBackward<Dimensions::Transform::TRA2D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2DTRA);

      // Compute projection transform for second dimension 
      coord.transform2D().project(rOutVar.rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA2D>::Type::ProjectorType::Id>(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2DTRA);

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(*pInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(*pInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

}
}
