/** 
 * @file BackwardSerialConfigurator.hpp
 * @brief This defines the backward transform serial operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BACKWARDSERIALCONFIGURATOR_HPP
#define BACKWARDSERIALCONFIGURATOR_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/TransformTreeSelector.hpp"
#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformConfigurators/BackwardConfigurator3D.hpp"
#else
   #include "TransformConfigurators/BackwardConfigurator2D.hpp"
#endif //GEOMHDISCC_SPATIALDIMENSION_3D

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the backward transform serial operations
    */
   class BackwardSerialConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::NONE;

         /**
          * @brief Compute the first step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void firstStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the second step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void secondStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Compute the last step in the backward transform
          *
          * @param tree       Transform projector tree
          * @param rVariable  Variable corresponding to the name
          * @param coord      Transform coordinator
          *
          * \tparam TVariable Type of the physical variable
          */
         template <typename TVariable> static void lastStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord);

         /**
          * @brief Setup first exchange communication
          */
         static void setup1DCommunication(const int packs, TransformCoordinatorType& coord);

         /**
          * @brief Setup second exchange communication
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
         BackwardSerialConfigurator();

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSerialConfigurator();

      private:
   };

   inline void BackwardSerialConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSerialConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSerialConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSerialConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSerialConfigurator::firstStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the transforms
      ProjectorSpecEdge_iterator itSpec;
      ProjectorPhysEdge_iterator itPhys;

      // Ranges for the vector of edges for the three transforms
      ProjectorSpecEdge_range rangeSpec = tree.edgeRange();
      ProjectorPhysEdge_range rangePhys;

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);

      #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
         // Iterators for the second transforms
         ProjectorPartEdge_iterator it2D;

         // Ranges for the vector of edges for the second transforms
         ProjectorPartEdge_range range2D;

         // Loop over first transform
         int holdSpec = std::distance(rangeSpec.first, rangeSpec.second) - 1;
         for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec, --holdSpec)
         {
            // Compute first transform
            BackwardConfigurator3D::project1D(*itSpec, coord, holdSpec);

            range2D = itSpec->edgeRange();
            int recover2D = 0;
            int hold2D = std::distance(range2D.first, range2D.second) - 1;
            for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
            {
               // Compute second transform
               BackwardConfigurator3D::project2D(*it2D, coord, recover2D, hold2D);

               rangePhys = it2D->edgeRange();
               int recoverPhys = 0;
               int holdPhys = std::distance(rangePhys.first, rangePhys.second) - 1;
               for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys, ++recoverPhys, --holdPhys)
               {
                  // Prepare physical output data
                  BackwardConfigurator3D::preparePhysical(tree, *itPhys, rVariable, coord);

                  // Compute third transform
                  BackwardConfigurator3D::projectND(*itPhys, coord, recoverPhys, holdPhys);
               }
            }
         }
      #else
         // Loop over first transform
         int holdSpec = std::distance(rangeSpec.first, rangeSpec.second) - 1;
         for(itSpec = rangeSpec.first; itSpec != rangeSpec.second; ++itSpec, --holdSpec)
         {
            // Compute first transform
            BackwardConfigurator2D::project1D(*itSpec, coord, holdSpec);

            rangePhys = itSpec->edgeRange();
            int recoverPhys = 0;
            int holdPhys = std::distance(rangePhys.first, rangePhys.second) - 1;
            for(itPhys = rangePhys.first; itPhys != rangePhys.second; ++itPhys, ++recoverPhys, --holdPhys)
            {
               // Prepare physical output data
               BackwardConfigurator2D::preparePhysical(tree, *itPhys, rVariable, coord);

               // Compute third transform
               BackwardConfigurator2D::projectND(*itPhys, coord, recoverPhys, holdPhys);
            }
         }
      #endif //GEOMHDISCC_SPATIALDIMENSION_3D
   }

   template <typename TVariable> void BackwardSerialConfigurator::secondStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSerialConfigurator::lastStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

}
}

#endif // BACKWARDSERIALCONFIGURATOR_HPP
