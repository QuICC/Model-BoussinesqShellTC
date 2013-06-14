/** \file ForwardSingle1DConfigurator.cpp
 *  \brief Source of the implementation of the forward single splitting algorithm on first dimension configurator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/ForwardSingle1DConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void ForwardSingle1DConfigurator::firstStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Only compute for equations requiring nonlinear terms
      if(spEquation->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear())
      {
         // Compute the nonlinear interaction
         ForwardConfigurator::nonlinearTerm<FieldComponents::Physical::SCALAR>(spEquation, coord);

         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute scalar integration in the third dimension
         ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_SCALAR>(coord);

         // Compute scalar integration in the second dimension
         ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_SCALAR>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
      }
   }

   void ForwardSingle1DConfigurator::firstStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Only compute for equations requiring nonlinear terms
      if(spEquation->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear())
      {
         // Compute the nonlinear interaction
         ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_ONE>(spEquation, coord);


         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in the third dimension
         ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_ONE>(coord);

         // Compute integration in the second dimension
         ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_ONE>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


         // Compute the nonlinear interaction
         ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_TWO>(spEquation, coord);

         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in the third dimension
         ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_TWO>(coord);

         // Compute integration in the second dimension
         ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_TWO>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


         // Compute the nonlinear interaction
         ForwardConfigurator::nonlinearTerm<TransformSteps::Physical::NONLINEAR_THREE>(spEquation, coord);

         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in third dimension
         ForwardConfigurator::integrate3D<TransformSteps::Forward<Dimensions::Transform::TRA3D>::STEP_VECTOR_THREE>(coord);

         // Compute integration in second dimension
         ForwardConfigurator::integrate2D<TransformSteps::Forward<Dimensions::Transform::TRA2D>::STEP_VECTOR_THREE>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
      }
   }

   void ForwardSingle1DConfigurator::lastStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Only compute for equations requiring nonlinear terms
      if(spEquation->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear())
      {
         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute scalar integration in the first dimension
         ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_SCALAR>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
      }
   }

   void ForwardSingle1DConfigurator::lastStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Only compute for equations requiring nonlinear terms
      if(spEquation->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear())
      {
         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in the second dimension for second component
         ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_ONE>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in the first dimension for third component
         ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_TWO>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


         // Start profiler
         ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

         // Compute integration in first dimension for first component
         ForwardConfigurator::integrate1D<TransformSteps::Forward<Dimensions::Transform::TRA1D>::STEP_VECTOR_THREE>(coord);

         // Stop profiler
         ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
      }
   }

   void ForwardSingle1DConfigurator::secondStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

   void ForwardSingle1DConfigurator::secondStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

}
}
