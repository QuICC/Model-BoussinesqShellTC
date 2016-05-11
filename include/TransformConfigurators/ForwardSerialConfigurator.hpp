/** 
 * @file ForwardSerialConfigurator.hpp
 * @brief This class defines the forward transform serial operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FORWARDSERIALCONFIGURATOR_HPP
#define FORWARDSERIALCONFIGURATOR_HPP

// Configuration includes
// 

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
    * @brief This class defines the forward transform serial operations
    */
   class ForwardSerialConfigurator: public ForwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::NONE;

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
         ForwardSerialConfigurator() {};

         /**
          * @brief Empty destructor
          */
         virtual ~ForwardSerialConfigurator() {};

      private:
   };

   inline void ForwardSerialConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
   }

   inline void ForwardSerialConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

   template <typename TSharedEquation> void ForwardSerialConfigurator::firstStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the transforms
      TransformTreeEdge::EdgeType_iterator itSpec;
      TransformTreeEdge::EdgeType_iterator itPhys;

      // Ranges for the vector of edges for the three transforms
      TransformTreeEdge::EdgeType_range rangeSpec;
      TransformTreeEdge::EdgeType_range rangePhys = tree.root().edgeRange();

      // Compute the nonlinear interaction
      ForwardConfigurator::nonlinearTerm(tree, spEquation, coord);

      // Start profiler
      ProfilerMacro_start(ProfilerMacro::FWDTRANSFORM);

      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         // Iterators for the second transforms
         TransformTreeEdge::EdgeType_iterator it2D;

         // Ranges for the vector of edges for the second transforms
         TransformTreeEdge::EdgeType_range range2D;

         // Loop over first transform
         int holdPhys = std::distance(rangePhys.first, rangePhys.second) - 1;
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys, --holdPhys)
         {
            // Compute third transform
            ForwardConfigurator3D::integrateND(*itPhys, coord, holdPhys);

            range2D = itPhys->edgeRange();
            int recover2D = 0;
            int hold2D = std::distance(range2D.first, range2D.second) - 1;
            for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
            {
               // Compute second transform
               ForwardConfigurator3D::integrate2D(*it2D, coord, recover2D, hold2D);

               rangeSpec = it2D->edgeRange();
               int recoverSpec = 0;
               int holdSpec = std::distance(rangeSpec.first, rangeSpec.second) - 1;
               for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec, ++recoverSpec, --holdSpec)
               {
                  // Compute third transform
                  ForwardConfigurator3D::integrate1D(*itSpec, coord, recoverSpec, holdSpec);

                  // Update equation
                  ForwardConfigurator3D::updateEquation(*itSpec, spEquation, coord, holdSpec);
               }
            }
         }
      #else
         // Loop over physical transform
         int holdPhys = std::distance(rangePhys.first, rangePhys.second) - 1;
         for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys, --holdPhys)
         {
            // Compute physical transform
            ForwardConfigurator2D::integrateND(*itPhys, coord, holdPhys);

            rangeSpec = itPhys->edgeRange();
            int recoverSpec = 0;
            int holdSpec = std::distance(rangeSpec.first, rangeSpec.second) - 1;
            for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec, ++recoverSpec, --holdSpec)
            {
               // Compute third transform
               ForwardConfigurator2D::integrate1D(*itSpec, coord, recoverSpec, holdSpec);

               // Update equation
               ForwardConfigurator2D::updateEquation(*itSpec, spEquation, coord, holdSpec);
            }
         }
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D

      // Stop profiler
      ProfilerMacro_stop(ProfilerMacro::FWDTRANSFORM);
   }

   template <typename TSharedEquation> void ForwardSerialConfigurator::secondStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }
   
   template <typename TSharedEquation> void ForwardSerialConfigurator::lastStep(const TransformTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a last step
   }

}
}

#endif // FORWARDSERIALCONFIGURATOR_HPP
