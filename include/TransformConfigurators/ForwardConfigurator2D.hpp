/** 
 * @file ForwardConfigurator2D.hpp
 * @brief This class defines the base operations for a forward transform in 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FORWARDCONFIGURATOR2D_HPP
#define FORWARDCONFIGURATOR2D_HPP

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
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/TransformStepsMacro.h"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the base operations for a forward transform in 2D space
    */
   class ForwardConfigurator2D
   {
      public:

      protected:
         /**
          * @brief Compute nonlinear interaction on a scalar or vector field
          *
          * @param spEquation Equation providing the nonlinear computation
          * @param coord      Transform coordinator
          */
         static void nonlinearTerm(const TransformTree& tree, Equations::SharedIEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the first dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrate1D(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Compute the integration transform of the last dimension
          *
          * @param coord   Transform coordinator
          */
         static void integrateND(const TransformTreeEdge& edge, TransformCoordinatorType& coord);

         /**
          * @brief Update variable values from dealiased data
          *
          * @param spEquation Equation providing the variable
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void updateEquation(const TransformTreeEdge& edge, TSharedEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Empty constructor
          */
         ForwardConfigurator2D() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardConfigurator2D() {};

      private: 
   };

   template <typename TSharedEquation> void ForwardConfigurator2D::updateEquation(const TransformTreeEdge& edge, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::DIAGNOSTICEQUATION);

      // Recover temporary storage
      TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute linear term component
      spEquation->updateDealiasedUnknown(rInVar, edge.outId<FieldComponents::Spectral::Id>(), edge.arithId());

      // Hold temporary storage
      if(edge.combinedOutId() >= 0)
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().holdBwd(rInVar, edge.combinedOutId());

      // Free the temporary storage
      } else
      {
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::DIAGNOSTICEQUATION);
   }

#ifdef QUICC_SPATIALDIMENSION_2D
   typedef ForwardConfigurator2D ForwardConfigurator;
#endif //QUICC_SPATIALDIMENSION_2D

}
}

#endif // FORWARDCONFIGURATOR2D_HPP
