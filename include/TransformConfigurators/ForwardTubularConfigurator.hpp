/** \file ForwardTubularConfigurator.hpp
 *  \brief This class defines the forward transform tubular splitting operations
 */
#ifdef GEOMHDISCC_MPIALGO_TUBULAR

#ifndef FORWARDTUBULARCONFIGURATOR_HPP
#define FORWARDTUBULARCONFIGURATOR_HPP

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
    * @brief This class defines the forward transform tubular splitting operations
    */
   class ForwardTubularConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::BOTH;

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
         ForwardTubularConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardTubularConfigurator() {};

      private:
   };

   inline void ForwardTubularConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs);

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardCommunication();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateForwardCommunication();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

}
}

#endif // FORWARDTUBULARCONFIGURATOR_HPP

#endif //GEOMHDISCC_MPIALGO_TUBULAR
