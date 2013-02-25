/** \file BackwardSingle1DConfigurator.cpp
 *  \brief Source of the implementation of the backward transform first single splitting algorith configurator
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "TransformConfigurators/BackwardSingle1DConfigurator.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace Transform {

   void BackwardSingle1DConfigurator::firstPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection(rScalar, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical(rScalar, coord);


      // Compute first backward transform
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::SCALAR>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection(rScalar, coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Phys::GRAD1D>(rScalar, coord);

      // Compute first backward transform for first gradient component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::GRAD1D>(coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Phys::GRAD2D>(rScalar, coord);

      // Compute first backward transform for second gradient component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::GRAD2D>(coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Phys::GRAD3D>(rScalar, coord);

      // Compute first backward transform for third gradient component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::GRAD3D>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECVECT1D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Phys::VECTOR1D>(rVector, coord);

      // Compute first backward transform for the first component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::VECTOR1D>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECVECT2D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Phys::VECTOR2D>(rVector, coord);

      // Compute first backward transform for the second component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::VECTOR2D>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECVECT3D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Phys::VECTOR3D>(rVector, coord);

      // Compute first backward transform for the third component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::VECTOR3D>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECCURL1D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Phys::CURL1D>(rVector, coord);

      // Compute first backward transform for first curl component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::CURL1D>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECCURL2D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Phys::CURL2D>(rVector, coord);

      // Compute first backward transform for second curl component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::CURL2D>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Bwd1D::SPECCURL3D>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Phys::CURL3D>(rVector, coord);

      // Compute first backward transform for third curl component
      BackwardConfigurator::project1D<TransformSteps::Bwd1D::CURL3D>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Compute second backward transform
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::SCALAR>(coord);

      // Compute third backward transform
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::SCALAR>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first gradient component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::GRAD1D>(coord);

      // Compute third backward transform for first gradient component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::GRAD1D>(coord);


      // Compute second backward transform for second gradient component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::GRAD2D>(coord);

      // Compute third backward transform of for second gradient component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::GRAD2D>(coord);


      // Compute second backward transform for third  gradient component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::GRAD3D>(coord);

      // Compute third backward transform for third gradient component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::GRAD3D>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::VECTOR1D>(coord);

      // Compute third backward transform for first component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::VECTOR1D>(coord);


      // Compute second backward transform for second component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::VECTOR2D>(coord);

      // Compute third backward transform for second component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::VECTOR2D>(coord);


      // Compute second backward transform for third component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::VECTOR3D>(coord);

      // Compute third backward transform for third component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::VECTOR3D>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first curl component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::CURL1D>(coord);

      // Compute third backward transform for first curl component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::CURL1D>(coord);


      // Compute second backward transform for second curl component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::CURL2D>(coord);

      // Compute third backward transform for second curl component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::CURL2D>(coord);


      // Compute second backward transform for third curl component
      BackwardConfigurator::project2D<TransformSteps::Bwd2D::CURL3D>(coord);

      // Compute third backward transform for third curl component
      BackwardConfigurator::project3D<TransformSteps::Bwd3D::CURL3D>(coord);
   }

   void BackwardSingle1DConfigurator::secondPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

   void BackwardSingle1DConfigurator::secondPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

   void BackwardSingle1DConfigurator::secondPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

   void BackwardSingle1DConfigurator::secondPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // No need for second step
   }

}
}
