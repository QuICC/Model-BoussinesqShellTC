/** 
 * @file ForwardTubularConfigurator.hpp
 * @brief This class defines the forward transform tubular splitting operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, true);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareForwardReceive();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, true);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().prepareForwardReceive();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateBackwardSend();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateBackwardSend();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardTubularConfigurator::firstStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorTree::Integrator3DEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorTree::Integrator3DEdge_range range3D = tree.edgeRange();

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
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardTubularConfigurator::secondStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorTree::Integrator2DEdge_iterator it2D;
      IntegratorTree::Integrator3DEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorTree::Integrator2DEdge_range range2D;
      IntegratorTree::Integrator3DEdge_range range3D = tree.edgeRange();

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
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
   
   template <typename TSharedEquation> void ForwardTubularConfigurator::lastStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorTree::Integrator1DEdge_iterator it1D;
      IntegratorTree::Integrator2DEdge_iterator it2D;
      IntegratorTree::Integrator3DEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorTree::Integrator1DEdge_range range1D;
      IntegratorTree::Integrator2DEdge_range range2D;
      IntegratorTree::Integrator3DEdge_range range3D = tree.edgeRange();

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

#endif // FORWARDTUBULARCONFIGURATOR_HPP

#endif //GEOMHDISCC_MPIALGO_TUBULAR
