/** 
 * @file BackwardSingle1DConfigurator.hpp
 * @brief This defines the backward transform first exchange single splitting operations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#if defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_COUPLED2D

#ifndef BACKWARDSINGLE1DCONFIGURATOR_HPP
#define BACKWARDSINGLE1DCONFIGURATOR_HPP

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
#include "TransformConfigurators/BackwardConfigurator3D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the backward transform first exchange single splitting operations
    */
   class BackwardSingle1DConfigurator: public BackwardConfigurator3D
   {
      public:
         /**
          * @brief Location of the splitting
          */
         static const Splitting::Locations::Id  SplitLocation = Splitting::Locations::FIRST;

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
         BackwardSingle1DConfigurator(){};

         /**
          * @brief Empty destructor
          */
         virtual ~BackwardSingle1DConfigurator(){};

      private:
   };

   template <typename TVariable> void BackwardSingle1DConfigurator::firstStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      ProjectorSpecEdge_iterator it1D;

      // Ranges for the vector of edges for the three transforms
      ProjectorSpecEdge_range range1D = tree.edgeRange();

      // Prepare required spectral data
      BackwardConfigurator3D::prepareSpectral(tree, rVariable, coord);

      // Loop over first transform
      int hold1D = std::distance(range1D.first, range1D.second) - 1;
      for(it1D = range1D.first; it1D != range1D.second; ++it1D, --hold1D)
      {
         // Compute first transform
         BackwardConfigurator3D::project1D(*it1D, coord, hold1D);
      }
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::secondStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
   }

   template <typename TVariable> void BackwardSingle1DConfigurator::lastStep(const ProjectorTree& tree, TVariable& rVariable, TransformCoordinatorType& coord)
   {
      // Iterators for the three transforms
      ProjectorSpecEdge_iterator it1D;
      ProjectorPartEdge_iterator it2D;
      ProjectorPhysEdge_iterator it3D;

      // Ranges for the vector of edges for the three transforms
      ProjectorSpecEdge_range range1D = tree.edgeRange();
      ProjectorPartEdge_range range2D;
      ProjectorPhysEdge_range range3D;

      // Loop over first transform
      for(it1D = range1D.first; it1D != range1D.second; ++it1D)
      {
         range2D = it1D->edgeRange();
         int recover2D = 0;
         int hold2D = std::distance(range2D.first, range2D.second) - 1;
         for(it2D = range2D.first; it2D != range2D.second; ++it2D, ++recover2D, --hold2D)
         {
            // Compute second transform
            BackwardConfigurator3D::project2D(*it2D, coord, recover2D, hold2D);

            range3D = it2D->edgeRange();
            int recover3D = 0;
            int hold3D = std::distance(range3D.first, range3D.second) - 1;
            for(it3D = range3D.first; it3D != range3D.second; ++it3D, ++recover3D, --hold3D)
            {
               // Prepare physical output data
               BackwardConfigurator3D::preparePhysical(tree, *it3D, rVariable, coord);

               // Compute third transform
               BackwardConfigurator3D::projectND(*it3D, coord, recover3D, hold3D);
            }
         }
      }
   }

   inline void BackwardSingle1DConfigurator::setup1DCommunication(const int packs, TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA2D>().setupCommunication(packs, TransformDirection::BACKWARD);

      coord.communicator().converter<Dimensions::Transform::TRA2D>().prepareBackwardReceive();
   }

   inline void BackwardSingle1DConfigurator::setup2DCommunication(const int packs, TransformCoordinatorType& coord)
   {
   }

   inline void BackwardSingle1DConfigurator::initiate1DCommunication(TransformCoordinatorType& coord)
   {
      coord.communicator().converter<Dimensions::Transform::TRA2D>().initiateForwardSend();
   }

   inline void BackwardSingle1DConfigurator::initiate2DCommunication(TransformCoordinatorType& coord)
   {
   }

}
}

#endif // BACKWARDSINGLE1DCONFIGURATOR_HPP

#endif //defined GEOMHDISCC_MPIALGO_SINGLE1D  || defined GEOMHDISCC_MPIALGO_COUPLED2D
