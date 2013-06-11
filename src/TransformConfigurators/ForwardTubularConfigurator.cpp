/** \file ForwardTubularConfigurator.cpp
 *  \brief Source of the implementation of the forward transform tubula algorithm configurator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/ForwardTubularConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardTubularConfigurator::firstStep(Equations::SharedIScalarPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<FieldComponents::Physical::SCALAR>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute scalar integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_SCALAR>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardTubularConfigurator::firstStep(Equations::SharedIVectorPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_ONE>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_ONE>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_TWO>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_TWO>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_THREE>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_THREE>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardTubularConfigurator::secondStep(Equations::SharedIScalarPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute scalar integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_SCALAR>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardTubularConfigurator::secondStep(Equations::SharedIVectorPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_ONE>(coord);


      // Compute integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_TWO>(coord);


      // Compute integration in second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_THREE>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardTubularConfigurator::lastStep(Equations::SharedIScalarPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute scalar integration in the first dimension
      ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_SCALAR>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardTubularConfigurator::lastStep(Equations::SharedIVectorPEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the second dimension
      ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_ONE>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the first dimension
      ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_TWO>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in first dimension
      ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_THREE>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

}
}
