/** 
 * @file BackwardSingle1DConfigurator.cpp
 * @brief Source of the implementation of the backward transform first single splitting algorith configurator
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
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_SCALAR>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection(rScalar, coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Physical::GRAD_ONE>(rScalar, coord);

      // Compute first backward transform for first gradient component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_GRAD_ONE>(coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Physical::GRAD_TWO>(rScalar, coord);

      // Compute first backward transform for second gradient component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_GRAD_TWO>(coord);


      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareGradient<TransformSteps::Physical::GRAD_THREE>(rScalar, coord);

      // Compute first backward transform for third gradient component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_GRAD_THREE>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECTOR_ONE>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Physical::VECTOR_ONE>(rVector, coord);

      // Compute first backward transform for the first component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_VECTOR_ONE>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECTOR_TWO>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Physical::VECTOR_TWO>(rVector, coord);

      // Compute first backward transform for the second component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_VECTOR_TWO>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECTOR_THREE>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::preparePhysical<TransformSteps::Physical::VECTOR_THREE>(rVector, coord);

      // Compute first backward transform for the third component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_VECTOR_THREE>(coord);
   }

   void BackwardSingle1DConfigurator::firstPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECURL_ONE>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Physical::CURL_ONE>(rVector, coord);

      // Compute first backward transform for first curl component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_CURL_ONE>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECURL_TWO>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Physical::CURL_TWO>(rVector, coord);

      // Compute first backward transform for second curl component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_CURL_TWO>(coord);


      // Prepare computation of projection
      BackwardConfigurator::prepareProjection<TransformSteps::Backward<Dimensions::Transform::TRA1D>::SPECURL_THREE>(rVector, coord);

      // Prepare computation of nonlinear interactions
      BackwardConfigurator::prepareCurl<TransformSteps::Physical::CURL_THREE>(rVector, coord);

      // Compute first backward transform for third curl component
      BackwardConfigurator::project1D<TransformSteps::Backward<Dimensions::Transform::TRA1D>::STEP_CURL_THREE>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Compute second backward transform
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_SCALAR>(coord);

      // Compute third backward transform
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_SCALAR>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysicalDiff(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first gradient component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_GRAD_ONE>(coord);

      // Compute third backward transform for first gradient component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_GRAD_ONE>(coord);


      // Compute second backward transform for second gradient component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_GRAD_TWO>(coord);

      // Compute third backward transform of for second gradient component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_GRAD_TWO>(coord);


      // Compute second backward transform for third  gradient component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_GRAD_THREE>(coord);

      // Compute third backward transform for third gradient component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_GRAD_THREE>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_VECTOR_ONE>(coord);

      // Compute third backward transform for first component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_VECTOR_ONE>(coord);


      // Compute second backward transform for second component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_VECTOR_TWO>(coord);

      // Compute third backward transform for second component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_VECTOR_TWO>(coord);


      // Compute second backward transform for third component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_VECTOR_THREE>(coord);

      // Compute third backward transform for third component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_VECTOR_THREE>(coord);
   }

   void BackwardSingle1DConfigurator::lastPhysicalDiff(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Compute second backward transform for first curl component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_CURL_ONE>(coord);

      // Compute third backward transform for first curl component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_CURL_ONE>(coord);


      // Compute second backward transform for second curl component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_CURL_TWO>(coord);

      // Compute third backward transform for second curl component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_CURL_TWO>(coord);


      // Compute second backward transform for third curl component
      BackwardConfigurator::project2D<TransformSteps::Backward<Dimensions::Transform::TRA2D>::STEP_CURL_THREE>(coord);

      // Compute third backward transform for third curl component
      BackwardConfigurator::project3D<TransformSteps::Backward<Dimensions::Transform::TRA3D>::STEP_CURL_THREE>(coord);
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
