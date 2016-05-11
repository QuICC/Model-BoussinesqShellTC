/** 
 * @file ForwardConfigurator2D.cpp
 * @brief Source of the implementation of the base forward configurator in 2D space
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
#include "TransformConfigurators/ForwardConfigurator2D.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardConfigurator2D::nonlinearTerm(const TransformTree& tree, Equations::SharedIEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::NONLINEAR);

      // Get physical storage
      TransformCoordinatorType::CommunicatorType::FwdNDType &rNLComp = coord.communicator().providePhysical();

      // Compute nonlinear term component
      spEquation->computeNonlinear(rNLComp, tree.comp<FieldComponents::Physical::Id>());
      spEquation->useNonlinear(rNLComp, tree.comp<FieldComponents::Physical::Id>());

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(rNLComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::NONLINEAR);
   }

   void ForwardConfigurator2D::integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate ND", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDND);

      // Get recover input data from hold
      TransformCoordinatorType::CommunicatorType::FwdNDType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRAND>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::BwdNDType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRAND>().provideBwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDNDTRA);

      // Compute integration transform of third dimension
      coord.transformND().integrate(rOutVar.rData(), rInVar.data(), edge.opId<TransformSelector<Dimensions::Transform::TRAND>::Type::IntegratorType::Id>(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDNDTRA);

      // Hold temporary input storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().holdFwd(rInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRAND>().freeFwd(rInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferBackward<Dimensions::Transform::TRAND>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDND);
   }

   void ForwardConfigurator2D::integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
      // Debugger message
      DebuggerMacro_msg("Integrate 1D", 4);

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1D);

      TransformCoordinatorType::CommunicatorType::Fwd1DType* pInVar;

      // Recover hold input data
      if(recover)
      {
         pInVar = &coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverFwd();

      // Get the transfered input data
      } else
      {
         pInVar = &coord.communicator().receiveForward<Dimensions::Transform::TRA1D>();
      }

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideBwd();

      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1DTRA);

      // Compute integration transform of first dimension
      coord.transform1D().integrate(rOutVar.rData(), pInVar->data(), edge.opId<TransformSelector<Dimensions::Transform::TRA1D>::Type::IntegratorType::Id>(), Arithmetics::SET);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1DTRA);

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdFwd(*pInVar);

      // Free temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(*pInVar);
      }

      // Hold temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1D);
   }

}
}
