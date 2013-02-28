/** \file BackwardConfigurator.cpp
 *  \brief Source of the implementation of the base backward configurator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/BackwardConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   template <> void BackwardConfigurator::prepareProjection<FieldComponents::Spectral::NOTUSED>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::prepareGradient<FieldComponents::Physical::NOTUSED>(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::preparePhysical<FieldComponents::Physical::NOTUSED>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::prepareCurl<FieldComponents::Physical::NOTUSED>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
   }

   void BackwardConfigurator::prepareProjection(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Put scalar into temporary hold storage
      coord.communicator().holdSpectral(rScalar.rDom(0).rTotal());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   void BackwardConfigurator::preparePhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Put scalar into temporary hold storage
      coord.communicator().holdPhysical(rScalar.rDom(0).rPhys());

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::Id::Bwd1D::NOTUSED>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage1D().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage1D().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.rData(), TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Free spectral input
      coord.communicator().freeSpectral(rInVar);

      // Transfer output data to next step
      coord.communicator().transferFwd1D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::DO_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage1D().recoverBwd();

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage1D().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.rData(), TransformCoordinatorType::Transform1DType::ProjectorType::DIFF);

      // Hold spectral input
      coord.communicator().holdSpectral(rInVar);

      // Transfer output data to next step
      coord.communicator().transferFwd1D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::Id::Bwd2D::NOTUSED>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rInVar = coord.communicator().receiveBwd2D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutVar = coord.communicator().storage2D().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.rData(), TransformCoordinatorType::Transform2DType::ProjectorType::PROJ);

      // Free temporary input storage
      coord.communicator().storage2D().freeBwd(rInVar);

      // Transfer output to next step
      coord.communicator().transferFwd2D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::DO_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd2DType &tpComp = coord.communicator().receiveBwd2D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rThComp = coord.communicator().storage2D().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rThComp.rData(), tpComp.rData(), TransformCoordinatorType::Transform2DType::ProjectorType::DIFF);

      // Transfer output data to next step
      coord.communicator().transferFwd2D(rThComp);

      // Hold temporary storage
      coord.communicator().storage2D().holdBwd(tpComp);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::Id::Bwd3D::NOTUSED>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rInVar = coord.communicator().receiveBwd3D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rOutVar = coord.communicator().storage3D().recoverFwd();

      // Compute projection transform for third dimension 
      coord.transform3D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.rData(), TransformCoordinatorType::Transform3DType::ProjectorType::PROJ);

      // Free temporary input storage
      coord.communicator().storage3D().freeBwd(rInVar);

      // Hold temporary storage
      coord.communicator().transferFwd3D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::BackwardBase::DO_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rInVar = coord.communicator().receiveBwd3D();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rOutVar = coord.communicator().storage3D().recoverFwd();

      // Compute projection transform for third dimension 
      coord.transform3D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.rData(), TransformCoordinatorType::Transform3DType::ProjectorType::DIFF);

      // Free temporary input storage
      coord.communicator().storage3D().freeBwd(rInVar);

      // Hold temporary storage
      coord.communicator().transferFwd3D(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

}
}
