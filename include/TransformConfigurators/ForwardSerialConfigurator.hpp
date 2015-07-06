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
#include "TypeSelectors/TransformTreeSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/ForwardConfigurator.hpp"

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

   template <typename TSharedEquation> void ForwardSerialConfigurator::firstStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      IntegratorSpecEdge_iterator it1D;
      IntegratorPartEdge_iterator it2D;
      IntegratorPhysEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      IntegratorSpecEdge_range range1D;
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

   template <typename TSharedEquation> void ForwardSerialConfigurator::secondStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a second step
   }
   
   template <typename TSharedEquation> void ForwardSerialConfigurator::lastStep(const IntegratorTree& tree, TSharedEquation spEquation, TransformCoordinatorType& coord)
   {
      // No need for a last step
   }

}
}

#endif // FORWARDSERIALCONFIGURATOR_HPP
