/** \file ForwardSingle1DConfigurator.hpp
 *  \brief This class defines the forward transform single first exchange splitting operations
 */
#if defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_FIXED

#ifndef FORWARDSINGLE1DCONFIGURATOR_HPP
#define FORWARDSINGLE1DCONFIGURATOR_HPP

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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward transform single first exchange splitting operations
    */
   class ForwardSingle1DConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::FIRST;

         /**
          * @brief First step in transform, including the nonlinear interaction for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void firstStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief First step in transform, including the nonlinear interaction for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void firstStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void secondStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void secondStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a scalar field
          *
          * @param spEquation Scalar equation
          * @param coord      Transform coordinator
          */
         static void lastStep(Equations::SharedIScalarEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform for a vector field
          *
          * @param spEquation Vector equation
          * @param coord      Transform coordinator
          */
         static void lastStep(Equations::SharedIVectorEquation spEquation, TransformCoordinatorType& coord);

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
         ForwardSingle1DConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSingle1DConfigurator() {};

      private:
   };

   inline void ForwardSingle1DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardSingle1DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSingle1DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardCommunication();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

}
}

#endif // FORWARDSINGLE1DCONFIGURATOR_HPP

#endif //defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_FIXED
