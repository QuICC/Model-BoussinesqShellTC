/** \file ForwardSingle2DConfigurator.hpp
    \brief This class defines the forward transform single second exchange splitting operations
 */

#ifndef FORWARDSINGLE2DCONFIGURATOR_HPP
#define FORWARDSINGLE2DCONFIGURATOR_HPP

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
#include "TransformConfigurators/ForwardConfigurator.hpp"

namespace EPMPhoenix {

   /**
    * @brief This class defines the forward transform single second exchange splitting operations
    */
   class ForwardSingle2DConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splittings::Locations::Id  SplitLocation = Splittings::Locations::SECOND;

         /**
          * @brief First step in transform, including the nonlinear interaction for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void firstStep(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief First step in transform, including the nonlinear interaction for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void firstStep(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void secondStep(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void secondStep(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void lastStep(SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void lastStep(SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief First exchange communication setup
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Second Exchange communication setup
          */
         static void setup2DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Initiate first exchange communication
          */
         static void initiate1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Initiate second exchange communication
          */
         static void initiate2DCommunication(TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Empty constructor
          */
         ForwardSingle2DConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSingle2DConfigurator() {};

      private:
   };

   inline void ForwardSingle2DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSingle2DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter3D().setupCommunication(packs);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardSingle2DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSingle2DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter3D().initiateForwardCommunication();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

}

#endif // FORWARDSINGLE2DCONFIGURATOR_HPP
