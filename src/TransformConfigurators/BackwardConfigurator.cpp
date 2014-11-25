/** 
 * @file BackwardConfigurator.cpp
 * @brief Source of the implementation of the base backward configurator
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
      coord.communicator().dealiasSpectral(rScalar.rDom(0).rTotal());

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

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Free spectral input
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::FINISH_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Free spectral input
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::START_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::ProjectorType::DIFF);

      // Hold spectral input
      coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

#if defined GEOMHDISCC_SPATIALSCHEME_CFT || GEOMHDISCC_SPATIALSCHEME_AFT || defined GEOMHDISCC_SPATIALSCHEME_WFT || defined GEOMHDISCC_SPATIALSCHEME_SLF || GEOMHDISCC_SPATIALSCHEME_BLF || defined GEOMHDISCC_SPATIALSCHEME_WLF

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::CONTINUE_DIVR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::ProjectorType::DIVR);

      // Hold spectral input
      coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <> void BackwardConfigurator::project1D<TransformSteps::BackwardBase::FINISH_DIVR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform1DType::ProjectorType::DIVR);

      // Free spectral input
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA1D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

#endif //defined GEOMHDISCC_SPATIALSCHEME_CFT || GEOMHDISCC_SPATIALSCHEME_AFT || defined GEOMHDISCC_SPATIALSCHEME_WFT || defined GEOMHDISCC_SPATIALSCHEME_SLF || GEOMHDISCC_SPATIALSCHEME_BLF || defined GEOMHDISCC_SPATIALSCHEME_WLF

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rInVar = coord.communicator().receiveBackward<Dimensions::Transform::TRA2D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform2DType::ProjectorType::PROJ);

      // Free temporary input storage
      coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(rInVar);

      // Transfer output to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::FINISH_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the input data from hold
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().recoverBwd();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform2DType::ProjectorType::PROJ);

      // Free temporary input storage
      coord.communicator().storage<Dimensions::Transform::TRA2D>().freeBwd(rInVar);

      // Transfer output to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA2D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   template <> void BackwardConfigurator::project2D<TransformSteps::BackwardBase::START_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD2D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd2DType &rInComp = coord.communicator().receiveBackward<Dimensions::Transform::TRA2D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd2DType &rOutComp = coord.communicator().storage<Dimensions::Transform::TRA2D>().provideFwd();

      // Compute projection transform for second dimension 
      coord.transform2D().project<Arithmetics::SET>(rOutComp.rData(), rInComp.data(), TransformCoordinatorType::Transform2DType::ProjectorType::DIFF);

      // Transfer output data to next step
      coord.communicator().transferForward<Dimensions::Transform::TRA2D>(rOutComp);

      // Hold temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA2D>().holdBwd(rInComp);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD2D);
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::BackwardBase::NOTHING>(TransformCoordinatorType& coord)
   {
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::BackwardBase::DO_SCALAR>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rInVar = coord.communicator().receiveBackward<Dimensions::Transform::TRA3D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA3D>().recoverFwd();

      // Compute projection transform for third dimension 
      coord.transform3D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform3DType::ProjectorType::PROJ);

      // Free temporary input storage
      coord.communicator().storage<Dimensions::Transform::TRA3D>().freeBwd(rInVar);

      // Hold temporary storage
      coord.communicator().transferForward<Dimensions::Transform::TRA3D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   template <> void BackwardConfigurator::project3D<TransformSteps::BackwardBase::DO_GRAD>(TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Get the transfered input data
      TransformCoordinatorType::CommunicatorType::Bwd3DType &rInVar = coord.communicator().receiveBackward<Dimensions::Transform::TRA3D>();

      // Get temporary storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA3D>().recoverFwd();

      // Compute projection transform for third dimension 
      coord.transform3D().project<Arithmetics::SET>(rOutVar.rData(), rInVar.data(), TransformCoordinatorType::Transform3DType::ProjectorType::DIFF);

      // Free temporary input storage
      coord.communicator().storage<Dimensions::Transform::TRA3D>().freeBwd(rInVar);

      // Hold temporary storage
      coord.communicator().transferForward<Dimensions::Transform::TRA3D>(rOutVar);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

}
}
