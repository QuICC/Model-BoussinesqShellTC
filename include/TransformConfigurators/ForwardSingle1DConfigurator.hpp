/** 
 * @file ForwardSingle1DConfigurator.hpp
 * @brief This class defines the forward transform single first exchange splitting operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#if defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_COUPLED2D

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
#include "TypeSelectors/TransformTreeSelector.hpp"
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
          * @brief First step in transform, including the nonlinear interaction
          *
          * @param spEquation Shared equation
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void firstStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform
          *
          * @param spEquation Shared equation
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void secondStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform
          *
          * @param spEquation Shared equation
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void lastStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

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

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareForwardReceive();

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

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateBackwardSend();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

   template <typename TSharedEquation> void ForwardSingle1DConfigurator::firstStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorPartEdge_iterator it2D;
      IntegratorPhysEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorPartEdge_range range2D;
      IntegratorPhysEdge_range range3D = tree.edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      int hold3D = std::distance(range3D.first, range3D.second) - 1;
      for(it3D = range3D.first; it3D != range3D.second; ++it3D, --hold3D)
      {
         // Compute third transform
         ForwardConfigurator::integrate3D(*it3D, coord, hold3D);

         range2D = it3D->edgeRange();
         int recover2D = 0;
         int hold2D = std::distance(range2D.first, range2D.second) - 1;
         for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
         {
            // Compute second transform
            ForwardConfigurator::integrate2D(*it2D, coord, recover2D, hold2D);
         }
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardSingle1DConfigurator::secondStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }
   
   template <typename TSharedEquation> void ForwardSingle1DConfigurator::lastStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorSpecEdge_iterator it1D;
      IntegratorPartEdge_iterator it2D;
      IntegratorPhysEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorSpecEdge_range range1D;
      IntegratorPartEdge_range range2D;
      IntegratorPhysEdge_range range3D = tree.edgeRange();

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         range2D = it3D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            range1D = it2D->edgeRange();
            int recover1D = 0;
            int hold1D = std::distance(range1D.first, range1D.second) - 1;
            for(it1D = range1D.first; it1D != range1D.second; ++it1D, ++recover1D, --hold1D)
            {
               // Compute third transform
               ForwardConfigurator::integrate1D(*it1D, coord, recover1D, hold1D);

               // Update equation
               ForwardConfigurator::updateEquation(*it1D, spEquation, coord, hold1D);
            }
         }
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

}
}

#endif // FORWARDSINGLE1DCONFIGURATOR_HPP

#endif //defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_COUPLED2D
