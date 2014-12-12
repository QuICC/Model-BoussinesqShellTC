/** 
 * @file ForwardConfigurator.cpp
 * @brief Source of the implementation of the base forward configurator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/ForwardConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardConfigurator::nonlinearTerm(const IntegratorTree& tree, Equations::SharedIEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::NONLINEAR);

      // Get physical storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rNLComp = coord.communicator().providePhysical();

      // Compute nonlinear term component
      spEquation->computeNonlinear(rNLComp, tree.comp());
      spEquation->useNonlinear(rNLComp, tree.comp());

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(rNLComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::NONLINEAR);
   }

   void ForwardConfigurator::integrate3D(const IntegratorTree::Integrator3DEdge& edge, TransformCoordinatorType& coord, const bool hold)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD3D);

      // Get recover input data from hold
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRA3D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA3D>().provideBwd();

      // Compute integration transform of third dimension
      coord.transform3D().integrate<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), edge.opId());

      // Hold temporary input storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA3D>().holdFwd(rInVar);

      // Free temporary input storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA3D>().freeFwd(rInVar);
      }

      // Transfer output data to next step
      coord.communicator().transferBackward<Dimensions::Transform::TRA3D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD3D);
   }

   void ForwardConfigurator::integrate2D(const IntegratorTree::Integrator2DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
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

      // Compute integration transform of second dimension
      coord.transform2D().integrate<Arithmetics::SET>(rOutVar.rData(), pInVar->data(), edge.opId());

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

   void ForwardConfigurator::integrate1D(const IntegratorTree::Integrator1DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold)
   {
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

      // Compute integration transform of first dimension
      coord.transform1D().integrate<Arithmetics::SET>(rOutVar.rData(), pInVar->data(), edge.opId());

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
