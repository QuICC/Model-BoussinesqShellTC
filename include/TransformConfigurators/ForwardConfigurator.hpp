/** \file ForwardConfigurator.hpp
 *  \brief This class defines the base operations for a forward transform
 *
 *  \mhdBug Needs test
 */

#ifndef FORWARDCONFIGURATOR_HPP
#define FORWARDCONFIGURATOR_HPP

// Configuration includes
// 
#include "Profiler/ProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformSteps.hpp"

namespace GeoMDHiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform
    */
   class ForwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Compute nonlinear interaction on a scalar
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          */
         static void nonlinearTerm(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute nonlinear interaction on a vector field
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          *
          * \tparam TComponent Physical vector field component
          */
         template <FieldComponents::Physical::Id TComponent> static void nonlinearTerm(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::::Id::Fwd1D::Step TStep> static void integrate1D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::::Id::Fwd2D::Step TStep> static void integrate2D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the third dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::::Id::Fwd3D::Step TStep> static void integrate3D(TransformCoordinatorType& coord);

         /**
          * @brief Compute linear interaction on a scalar
          *
          * @param spEquation Equation providing the linear computation
          * @param coord      Transform coordinator
          */
         static void linearTerm(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute linear interaction on a vector field
          *
          * @param spEquation Equation providing the linear computation
          * @param coord      Transform coordinator
          *
          * \tparam TComponent Spectral vector field component
          */
         template <FieldComponents::Spectral::Id TComponent> static void linearTerm(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Prepare the timestep RHS for a scalar
          *
          * @param spEquation Equation providing the timestep structure
          * @param coord      Transform coordinator
          */
         static void prepareTimestep(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Prepare the timestep RHS for a vector field
          *
          * @param spEquation Equation providing the timestep structure
          * @param coord      Transform coordinator
          *
          * \tparam TComponent Spectral vector field component
          */
         template <FieldComponents::Spectral::Id TComponent> static void prepareTimestep(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         ForwardConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardConfigurator() {};

      private: 
         template <Dimensions::Type, Dimensions::Type> TransformCoordinatorType::CommunicatorType::Fwd1DType& getFwd1D(TransformCoordinatorType& coord);
   };

   template <FieldComponents::Physical::Id TComponent> void ForwardConfigurator::nonlinearTerm(SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::NONLINEAR);

      // Get physical storage
      TransformCoordinatorType::CommunicatorType::Fwd3DType &rNLComp = coord.communicator().providePhysical();

      // Compute nonlinear term component
      spEquation->computeNonlinear(rNLComp, TComponent);

      // Transfer physical storage to next step
      coord.communicator().holdPhysical(rNLComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::NONLINEAR);
   }

   template <FieldComponents::Spectral::Id TComponent> void ForwardConfigurator::linearTerm(SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::LINEAR);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rComp = coord.communicator().storage1D().recoverBwd();

      // Compute linear term component
      spEquation->computeLinear(rComp, TComponent);

      // Hold the temporary storage
      coord.communicator().storage1D().holdBwd(rComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::LINEAR);
   }

   template <FieldComponents::Spectral::Id TComponent> void ForwardConfigurator::prepareTimestep(SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rComp = coord.communicator().storage1D().recoverBwd();

      // Compute linear term component
      spEquation->prepareTimestep(rComp, TComponent);

      // Free the temporary storage
      coord.communicator().storage1D().freeBwd(rComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);
   }

   /// Specialised 1D integration for NONE
   template <> void ForwardConfigurator::integrate1D<TransformSteps::::Id::Fwd1D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 1D integration for SCALAR
   template <> void ForwardConfigurator::integrate1D<TransformSteps::::Id::Fwd1D::SCALAR>(TransformCoordinatorType& coord);


   /// Specialised 2D integration for NONE
   template <> void ForwardConfigurator::integrate2D<TransformSteps::::Id::Fwd2D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 2D integration for SCALAR
   template <> void ForwardConfigurator::integrate2D<TransformSteps::::Id::Fwd2D::SCALAR>(TransformCoordinatorType& coord);


   /// Specialised 3D integration for NONE
   template <> void ForwardConfigurator::integrate3D<TransformSteps::::Id::Fwd3D::NONE>(TransformCoordinatorType& coord);

   /// Specialised 3D integration for SCALAR
   template <> void ForwardConfigurator::integrate3D<TransformSteps::::Id::Fwd3D::SCALAR>(TransformCoordinatorType& coord);


   /// Specialised linear term for NONE
   template <> void ForwardConfigurator::linearTerm<FieldComponents::Spectral::NONE>(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

   /// Specialised timestep preparation for NONE
   template <> void ForwardConfigurator::prepareTimestep<FieldComponents::Spectral::NONE>(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

}
}

#endif // FORWARDCONFIGURATOR_HPP
