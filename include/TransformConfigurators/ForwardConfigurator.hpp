/** \file ForwardConfigurator.hpp
 *  \brief This class defines the base operations for a forward transform
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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform
    */
   class ForwardConfigurator
   {
      public:
         /**
          * @brief Prepare the timestep RHS for a scalar
          *
          * @param spEquation Equation providing the timestep structure
          * @param coord      Transform coordinator
          */
         static void prepareTimestep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Prepare the timestep RHS for a vector
          *
          * @param spEquation Equation providing the timestep structure
          * @param coord      Transform coordinator
          */
         static void prepareTimestep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Compute nonlinear interaction on a scalar
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          */
         static void nonlinearTerm(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute nonlinear interaction on a vector field
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          *
          * \tparam TComponent Physical vector field component
          */
         template <FieldComponents::Physical::Id TComponent> static void nonlinearTerm(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::ForwardBase::Step TStep> static void integrate1D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::ForwardBase::Step TStep> static void integrate2D(TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the third dimension
          *
          * @param coord   Transform coordinator
          */
         template <TransformSteps::ForwardBase::Step TStep> static void integrate3D(TransformCoordinatorType& coord);

         /**
          * @brief Prepare the timestep RHS for a vector field
          *
          * @param spEquation Equation providing the timestep structure
          * @param coord      Transform coordinator
          *
          * \tparam TComponent Spectral vector field component
          */
         template <FieldComponents::Spectral::Id TComponent> static void prepareTimestep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         ForwardConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardConfigurator() {};

      private: 
   };

   template <FieldComponents::Physical::Id TComponent> void ForwardConfigurator::nonlinearTerm(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
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

   template <FieldComponents::Spectral::Id TComponent> void ForwardConfigurator::prepareTimestep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::TIMESTEP);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rComp = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute linear term component
      spEquation->prepareTimestep(rComp, TComponent);

      // Free the temporary storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rComp);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::TIMESTEP);
   }

   /// Specialised integration to do nothing
   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord);

   /// Specialised integration to compute scalar
   template <> void ForwardConfigurator::integrate1D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord);


   /// Specialised 2D integration to do nothing
   template <> void ForwardConfigurator::integrate2D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord);

   /// Specialised 2D integration to compute scalar
   template <> void ForwardConfigurator::integrate2D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord);


   /// Specialised 3D integration to do nothing
   template <> void ForwardConfigurator::integrate3D<TransformSteps::ForwardBase::NOTHING>(TransformCoordinatorType& coord);

   /// Specialised 3D integration to compute scalar
   template <> void ForwardConfigurator::integrate3D<TransformSteps::ForwardBase::DO_SCALAR>(TransformCoordinatorType& coord);


   /// Specialised timestep preparation to do nothing
   template <> void ForwardConfigurator::prepareTimestep<FieldComponents::Spectral::NOTUSED>(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

}
}

#endif // FORWARDCONFIGURATOR_HPP
