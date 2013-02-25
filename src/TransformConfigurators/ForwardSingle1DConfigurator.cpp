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

   void ForwardSingle1DConfigurator::firstStep(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute scalar integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Fwd3D::SCALAR>(coord);

      // Compute scalar integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Fwd2D::SCALAR>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardSingle1DConfigurator::firstStep(SharedVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Phys::NONLINEAR1D>(spEquation, coord);


      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Fwd3D::VECTOR1D>(coord);

      // Compute integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Fwd2D::VECTOR1D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Phys::NONLINEAR2D>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Fwd3D::VECTOR2D>(coord);

      // Compute integration in the second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Fwd2D::VECTOR2D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);


      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm<TransformSteps::Phys::NONLINEAR3D>(spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in third dimension
      ForwardConfigurator::integrate3D<TransformSteps::Fwd3D::VECTOR3D>(coord);

      // Compute integration in second dimension
      ForwardConfigurator::integrate2D<TransformSteps::Fwd2D::VECTOR3D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   void ForwardSingle1DConfigurator::lastStep(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute scalar integration in the first dimension
      ForwardConfigurator::integrate1D<TransformSteps::Fwd1D::SCALAR>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);

      // Compute the linear term
      ForwardConfigurator::linearTerm(spEquation, coord);

      // Prepare the timestep RHS
      ForwardConfigurator::prepareTimestep(spEquation, coord);
   }

   void ForwardSingle1DConfigurator::lastStep(SharedVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the second dimension for second component
      ForwardConfigurator::integrate1D<TransformSteps::Fwd1D::VECTOR1D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);

      // Add the toroidal linear term
      ForwardConfigurator::linearTerm<TransformSteps::Fwd1D::SPECVECT1D>(spEquation, coord);

      // Prepare the toroidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Fwd1D::SPECVECT1D>(spEquation, coord);


      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in the first dimension for third component
      ForwardConfigurator::integrate1D<TransformSteps::Fwd1D::VECTOR2D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);

      // Add the toroidal linear term
      ForwardConfigurator::linearTerm<TransformSteps::Fwd1D::SPECVECT2D>(spEquation, coord);

      // Prepare the toroidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Fwd1D::SPECVECT2D>(spEquation, coord);


      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Compute integration in first dimension for first component
      ForwardConfigurator::integrate1D<TransformSteps::Fwd1D::VECTOR3D>(coord);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);

      // Add the poloidal linear term
      ForwardConfigurator::linearTerm<TransformSteps::Fwd1D::SPECVECT3D>(spEquation, coord);

      // Prepare the poloidal timestep RHS
      ForwardConfigurator::prepareTimestep<TransformSteps::Fwd1D::SPECVECT3D>(spEquation, coord);
   }

   void ForwardSingle1DConfigurator::secondStep(SharedScalarEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

   void ForwardSingle1DConfigurator::secondStep(SharedVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

}
}
