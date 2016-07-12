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
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformConfigurators/ForwardConfigurator3D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward transform tubular splitting operations
    */
   class ForwardTubularConfigurator: public ForwardConfigurator3D
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

      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::FORWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareForwardReceive();

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   inline void ForwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, TransformDirection::FORWARD);

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

   template <typename TSharedEquation> void ForwardTubularConfigurator::firstStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator3D::nonlinearTerm(tree, spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         // Compute third transform
         ForwardConfigurator3D::integrateND(*it3D, coord);
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardTubularConfigurator::secondStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it2D;
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range2D;
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         range2D = it3D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            // Compute second transform
            ForwardConfigurator3D::integrate2D(*it2D, coord);
         }
      }

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }
   
   template <typename TSharedEquation> void ForwardTubularConfigurator::lastStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      TransformTreeEdge::EdgeType_citerator it1D;
      TransformTreeEdge::EdgeType_citerator it2D;
      TransformTreeEdge::EdgeType_citerator it3D;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_crange range1D;
      TransformTreeEdge::EdgeType_crange range2D;
      TransformTreeEdge::EdgeType_crange range3D = tree.root().edgeRange();

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      // Loop over first transform
      for(it3D = range3D.first; it3D != range3D.second; ++it3D)
      {
         range2D = it3D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            range1D = it2D->edgeRange();
            for(it1D = range1D.first; it1D != range1D.second; ++it1D)
            {
               // Compute third transform
               ForwardConfigurator3D::integrate1D(*it1D, coord);

               // Update equation
               ForwardConfigurator3D::updateEquation(*it1D, spEquation, coord);
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
