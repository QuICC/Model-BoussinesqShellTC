/** \file BackwardConfigurator.hpp
 *  \brief This class defines the base operations for a backward transform
 */

#ifndef BACKWARDCONFIGURATOR_HPP
#define BACKWARDCONFIGURATOR_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformSteps.hpp"

namespace GeoMDHiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a backward transform
    */
   class BackwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Prepare computation of projection for a scalar variable
          *
          * @param rScalar  Scalar variable
          * @param coord      Transform coordinator
          */
         static void prepareProjection(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of projection for a vector variable
          *
          * @param rVector  Vector variable
          * @param coord      Transform coordinator
          */
         template <FieldComponents::Spectral::Id TComponent> static void prepareProjection(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a scalar variable
          *
          * @param rScalar  Scalar variable
          * @param coord      Transform coordinator
          */
         static void preparePhysical(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of gradient projection for a scalar variable
          *
          * @param rScalar  Scalar variable
          * @param coord      Transform coordinator
          */
         template <FieldComponents::Physical::Id TComponent> static void prepareGradient(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of physical projection for a vector variable
          *
          * @param rVector  Vector variable
          * @param coord      Transform coordinator
          */
         template <FieldComponents::Physical::Id TComponent> static void preparePhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Prepare computation of curl projection for a vector variable
          *
          * @param rVector  Vector variable
          * @param coord      Transform coordinator
          */
         template <FieldComponents::Physical::Id TComponent> static void prepareCurl(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::Backward<Dimensions::Transform::TRA1D>::Step TStep> static void project1D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::Id::Bwd2D::Step TStep> static void project2D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the projection transform of the third dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::Id::Bwd3D::Step TStep> static void project3D(TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         BackwardConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardConfigurator() {};

      private: 
   };

   template <FieldComponents::Spectral::Id TComponent> void BackwardConfigurator::prepareProjection(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD1D);

      // Put scalar into temporary hold storage
      coord.communicator().holdSpectral(rVector.rDom(0).rTotal().rComp(static_cast<int>(TComponent)));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD1D);
   }

   template <FieldComponents::Physical::Id TComponent> void BackwardConfigurator::prepareGradient(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Put scalar into temporary hold storage
      coord.communicator().holdPhysical(rScalar.rDom(0).rGrad().rComp(static_cast<int>(TComponent)));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   template <FieldComponents::Physical::Id TComponent> void BackwardConfigurator::preparePhysical(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Put scalar into temporary hold storage
      coord.communicator().holdPhysical(rVector.rDom(0).rPhys().rComp(static_cast<int>(TComponent)));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   template <FieldComponents::Physical::Id TComponent> void BackwardConfigurator::prepareCurl(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWD3D);

      // Put scalar into temporary hold storage
      coord.communicator().holdPhysical(rVector.rDom(0).rCurl().rComp(static_cast<int>(TComponent)));

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWD3D);
   }

   /// Specialised projection preparation for NONE
   template <> void BackwardConfigurator::prepareProjection<FieldComponents::Spectral::NONE>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);


   /// Specialised gradient preparation for NONE
   template <> void BackwardConfigurator::prepareGradient<FieldComponents::Physical::NONE>(Datatypes::ScalarVariableType& rScalar, TransformCoordinatorType& coord);

   /// Specialised physical preparation for NONE
   template <> void BackwardConfigurator::preparePhysical<FieldComponents::Physical::NONE>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);

   /// Specialised curl preparation for NONE
   template <> void BackwardConfigurator::prepareCurl<FieldComponents::Physical::NONE>(Datatypes::VectorVariableType& rVector, TransformCoordinatorType& coord);


   /// Specialised 1D projection for NONE
   template <> void BackwardConfigurator::project1D<TransformSteps::Id::Bwd1D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 1D projection for SCALAR
   template <> void BackwardConfigurator::project1D<TransformSteps::Id::Bwd1D::SCALAR>(TransformCoordinatorType& coord);

   /// Specialised 1D projection for GRADX
   template <> void BackwardConfigurator::project1D<TransformSteps::Id::Bwd1D::GRADX>(TransformCoordinatorType& coord);


   /// Specialised 2D projection for NONE
   template <> void BackwardConfigurator::project2D<TransformSteps::Id::Bwd2D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 2D projection for SCALAR
   template <> void BackwardConfigurator::project2D<TransformSteps::Id::Bwd2D::SCALAR>(TransformCoordinatorType& coord);

   /// Specialised 2D projection for GRADY
   template <> void BackwardConfigurator::project2D<TransformSteps::Id::Bwd2D::GRADY>(TransformCoordinatorType& coord);


   /// Specialised 3D projection for NONE
   template <> void BackwardConfigurator::project3D<TransformSteps::Id::Bwd3D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 3D projection for SCALAR
   template <> void BackwardConfigurator::project3D<TransformSteps::Id::Bwd3D::SCALAR>(TransformCoordinatorType& coord);

   /// Specialised 3D projection for GRADZ
   template <> void BackwardConfigurator::project3D<TransformSteps::Id::Bwd3D::GRADZ>(TransformCoordinatorType& coord);

}
}

#endif // BACKWARDCONFIGURATOR_HPP
