/** \file ForwardConfigurator.cpp
 *  \brief Source of the implementation of the base forward configurator
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

   void ForwardConfigurator::prepareTimestep(Equations::SharedIScalarPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rScalar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute linear term component
      spEquation->prepareTimestep(rScalar, FieldComponents::Spectral::SCALAR);

      // Free the temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rScalar);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);
   }

   void ForwardConfigurator::prepareTimestep(Equations::SharedIVectorPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Prepare the toroidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Forward<Dimensions::Transform::TRA1D>::SPECTOR_ONE>(spEquation, coord);

      // Prepare the toroidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Forward<Dimensions::Transform::TRA1D>::SPECTOR_TWO>(spEquation, coord);

      // Prepare the poloidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Forward<Dimensions::Transform::TRA1D>::SPECTOR_THREE>(spEquation, coord);
   }

   template <> void ForwardConfigurator::prepareTimestep<FieldComponents::Spectral::NOTUSED>(Equations::SharedIVectorPEquation spEquation, TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRA1D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideBwd();

      // Compute integration transform of first dimension
      coord.transform1D().integrate<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Free temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rInVar);

      // Hold temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD1D);
   }

   template <> void ForwardConfigurator::integrate2D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::integrate2D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD2D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRA2D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideBwd();

      // Compute integration transform of second dimension
      coord.transform2D().integrate<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform2DType::IntegratorType::INTG);

      // Free temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA2D>().freeFwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferBackward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD2D);
   }

   template <> void ForwardConfigurator::integrate3D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::integrate3D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD3D);

      // Get recover input data from hold
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rInVar = coord.communicator().receiveForward<Dimensions::Transform::TRA3D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA3D>().provideBwd();

      // Compute integration transform of third dimension
      coord.transform3D().integrate<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform3DType::IntegratorType::INTG);

      // Free temporary input storage
      coord.communicator().storage<Dimensions::Transform::TRA3D>().freeFwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferBackward<Dimensions::Transform::TRA3D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD3D);
   }

}
}
