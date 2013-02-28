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

   template <> void ForwardConfigurator::linearTerm<FieldComponents::Spectral::NOTUSED>(SharedVectorEquation spEquation, TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::prepareTimestep<FieldComponents::Spectral::NOTUSED>(SharedVectorEquation spEquation, TransformCoordinatorType& coord)
   {
   }

   void ForwardConfigurator::nonlinearTerm(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::NONLINEAR);

      // Get physical storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rNLComp = coord.communicator().providePhysical();

      // Compute nonlinear term component
      spEquation->computeNonlinear(rNLComp);

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(rNLComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::NONLINEAR);
   }

   void ForwardConfigurator::linearTerm(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::LINEAR);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rScalar = coord.communicator().storage1D().recoverBwd();

      // Compute linear term component
      spEquation->computeLinear(rScalar);

      // Hold temporary storage
      coord.communicator().storage1D().holdBwd(rScalar);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::LINEAR);
   }

   void ForwardConfigurator::prepareTimestep(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rScalar = coord.communicator().storage1D().recoverBwd();

      // Compute linear term component
      spEquation->prepareTimestep(rScalar);

      // Free the temporary storage
      coord.communicator().storage1D().freeBwd(rScalar);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);
   }

   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWD1D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rInVar = coord.communicator().receiveFwd1D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rScalar = coord.communicator().storage1D().provideBwd();

      // Compute integration transform of first dimension
      coord.transform1D().integrate<Arithmetics::SET>(rScalar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

      // Free temporary storage
      coord.communicator().storage1D().freeFwd(rInVar);

      // Hold temporary storage
      coord.communicator().storage1D().holdBwd(rScalar);

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
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rInVar = coord.communicator().receiveFwd2D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rScalar = coord.communicator().storage2D().provideBwd();

      // Compute integration transform of second dimension
      coord.transform2D().integrate<Arithmetics::SET>(rScalar.rData(), rInVar.data(), TransformCoordinatorType::Transform2DType::IntegratorType::INTG);

      // Free temporary storage
      coord.communicator().storage2D().freeFwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferBwd2D(rScalar);

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
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rInVar = coord.communicator().receiveFwd3D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rOutVar = coord.communicator().storage3D().provideBwd();

      // Compute integration transform of third dimension
      coord.transform3D().integrate<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform3DType::IntegratorType::INTG);

      // Free temporary input storage
      coord.communicator().storage3D().freeFwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferBwd3D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWD3D);
   }

}
}
