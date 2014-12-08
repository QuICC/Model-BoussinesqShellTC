/** 
 * @file ForwardConfigurator.hpp
 * @brief This class defines the base operations for a forward transform
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "TransformConfigurators/TransformStepsMacro.h"
#include "TransformConfigurators/IntegratorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform
    */
   class ForwardConfigurator
   {
      public:

      protected:
         /**
          * @brief Compute nonlinear interaction on a scalar or vector field
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          */
         static void nonlinearTerm(const IntegratorTree& tree, Equations::SharedIEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate1D(const IntegratorTree::Integrator1DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold);

         /**
          * @brief Compute the integration transform of the second dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate2D(const IntegratorTree::Integrator2DEdge& edge, TransformCoordinatorType& coord, const bool recover, const bool hold);

         /**
          * @brief Compute the integration transform of the third dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate3D(const IntegratorTree::Integrator3DEdge& edge, TransformCoordinatorType& coord, const bool hold);

         /**
          * @brief Update variable values from dealiased data
          *
          * @param spEquation Equation providing the variable
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void updateEquation(const IntegratorTree::Integrator1DEdge& edge, TSharedEquation spEquation, TransformCoordinatorType& coord, const bool hold);

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

   template <typename TSharedEquation> void ForwardConfigurator::updateEquation(const IntegratorTree::Integrator1DEdge& edge, TSharedEquation spEquation, TransformCoordinatorType& coord, const bool hold)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute linear term component
      spEquation->updateDealiasedUnknown(rInVar, edge.specId());

      // Hold temporary storage
      if(hold)
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdFwd(rInVar);

      // Free the temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
   }

}
}

#endif // FORWARDCONFIGURATOR_HPP
