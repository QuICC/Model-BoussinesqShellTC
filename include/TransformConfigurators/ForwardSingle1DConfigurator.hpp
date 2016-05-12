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
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/ForwardConfiguratorMacro.h"

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
         template <typename TSharedEquation> static void firstStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Second step in transform
          *
          * @param spEquation Shared equation
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void secondStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

         /**
          * @brief Last step in transform
          *
          * @param spEquation Shared equation
          * @param coord      Transform coordinator
          */
         template <typename TSharedEquation> static void lastStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord);

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

   template <typename TSharedEquation> void ForwardSingle1DConfigurator::firstStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;

      // Ranges for the vector of edges for the transforms
      TransformTreeEdge::EdgeType_crange rangePhys = tree.root().edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         // Iterators for the second transforms
         IntegratorPartEdge_iterator it2D;

         // Ranges for the vector of edges for the second transforms
         IntegratorPartEdge_range range2D;

         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            // Compute third transform
            ForwardConfigurator3D::integrateND(*itPhys, coord);

            range2D = itPhys->edgeRange();
            for(it2D = range2D.first; it2D != range2D.second; ++it2D)
            {
               // Compute second transform
               ForwardConfigurator3D::integrate2D(*it2D, coord);
            }
         }
      #else
         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            // Compute third transform
            ForwardConfigurator2D::integrateND(*itPhys, coord);
         }
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardSingle1DConfigurator::secondStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }
   
   template <typename TSharedEquation> void ForwardSingle1DConfigurator::lastStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the transforms
      TransformTreeEdge::EdgeType_citerator itSpec;
      TransformTreeEdge::EdgeType_citerator itPhys;

      // Ranges for the vector of edges for the transforms
      TransformTreeEdge::EdgeType_crange rangeSpec;
      TransformTreeEdge::EdgeType_crange rangePhys = tree.root().edgeRange();

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         // Iterators for the second transforms
         IntegratorPartEdge_iterator it2D;

         // Ranges for the vector of edges for the second transforms
         IntegratorPartEdge_range range2D;

         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            range2D = itPhys->edgeRange();
            for(it2D = range2D.first; it2D != range2D.second; ++it2D)
            {
               rangeSpec = it2D->edgeRange();
               for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
               {
                  // Compute third transform
                  ForwardConfigurator3D::integrate1D(*itSpec, coord);

                  // Update equation
                  ForwardConfigurator3D::updateEquation(*itSpec, spEquation, coord);
               }
            }
         }
      #else
         // Loop over first transform
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys)
         {
            rangeSpec = itPhys->edgeRange();
            for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec)
            {
               // Compute third transform
               ForwardConfigurator2D::integrate1D(*itSpec, coord, recoverSpec);

               // Update equation
               ForwardConfigurator2D::updateEquation(*itSpec, spEquation, coord);
            }
         }
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

}
}

#endif // FORWARDSINGLE1DCONFIGURATOR_HPP

#endif //defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_COUPLED2D
