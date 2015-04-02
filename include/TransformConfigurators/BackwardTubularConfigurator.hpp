/** 
 * @file BackwardTubularConfigurator.hpp
 * @brief This defines the backward transform tubular splitting operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_MPIALGO_TUBULAR

#ifndef BACKWARDTUBULARCONFIGURATOR_HPP
#define BACKWARDTUBULARCONFIGURATOR_HPP

// Configuration includes
// 

// System includes
//

// External includes
//

// Project includes
//
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TransformConfigurators/BackwardConfigurator.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the backward transform tubular splitting operations
    */
   class BackwardTubularConfigurator: public BackwardConfigurator
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::BOTH;

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
         BackwardTubularConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardTubularConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardTubularConfigurator::firstStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      ProjectorTree::Projector1DEdge_iterator it1D;

      // Ranges for the vector of edges for the three transforms
      ProjectorTree::Projector1DEdge_range range1D = tree.edgeRange();

      // Prepare required spectral data
      BackwardConfigurator::prepareSpectral(tree, rVariable, coord);

      // Loop over first transform
      int hold1D = std::distance(range1D.first, range1D.second) - 1;
      for(it1D = range1D.first; it1D != range1D.second; ++it1D, --hold1D)
      {
         // Compute first transform
         BackwardConfigurator::project1D(*it1D, coord, hold1D);
      }
   }

   template <typename TVariable> void BackwardTubularConfigurator::secondStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      ProjectorTree::Projector1DEdge_iterator it1D;
      ProjectorTree::Projector2DEdge_iterator it2D;

      // Ranges for the vector of edges for the three transforms
      ProjectorTree::Projector1DEdge_range range1D = tree.edgeRange();
      ProjectorTree::Projector2DEdge_range range2D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         range2D = it1D->edgeRange();
         int recover2D = 0;
         int hold2D = std::distance(range2D.first, range2D.second) - 1;
         for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
         {
            // Compute second transform
            BackwardConfigurator::project2D(*it2D, coord, recover2D, hold2D);
         }
      }
   }

   template <typename TVariable> void BackwardTubularConfigurator::lastStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      ProjectorTree::Projector1DEdge_iterator it1D;
      ProjectorTree::Projector2DEdge_iterator it2D;
      ProjectorTree::Projector3DEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      ProjectorTree::Projector1DEdge_range range1D = tree.edgeRange();
      ProjectorTree::Projector2DEdge_range range2D;
      ProjectorTree::Projector3DEdge_range range3D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         range2D = it1D->edgeRange();
         for(it2D = range2D.first; it2D != range2D.second; ++it2D)
         {
            range3D = it2D->edgeRange();
            int recover3D = 0;
            int hold3D = std::distance(range3D.first, range3D.second) - 1;
            for(it3D = range3D.first; it3D != range3D.second; ++it3D, ++recover3D, --hold3D)
            {
               // Prepare physical output data
               BackwardConfigurator::preparePhysical(tree, *it3D, rVariable, coord);

               // Compute third transform
               BackwardConfigurator::project3D(*it3D, coord, recover3D, hold3D);
            }
         }
      }
   }

   inline void BackwardTubularConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, false);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareBackwardReceive();
   }

   inline void BackwardTubularConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA3D>().setupCommunication(packs, false);

      coord.communicator().converter<Dimensions::Transform::TRA3D>().prepareBackwardReceive();
   }

   inline void BackwardTubularConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardSend();
   }

   inline void BackwardTubularConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA3D>().initiateForwardSend();
   }

}
}

#endif // BACKWARDTUBULARCONFIGURATOR_HPP

#endif //GEOMHDISCC_MPIALGO_TUBULAR
