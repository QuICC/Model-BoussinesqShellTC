/** 
 * @file BackwardSingle2DGrouper.hpp
 * @brief This class defines the backward single grouping exchange grouping algorithm for the second exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef QUICC_TRANSGROUPER_SINGLE2D

#ifndef BACKWARDSINGLE2DGROUPER_HPP
#define BACKWARDSINGLE2DGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "TransformGroupers/IBackwardGrouper3D.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class BackwardSingle2DGrouper: public IBackwardGrouper3D
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardSingle2DGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardSingle2DGrouper();

         /**
          * @brief Setup the full backward transform structure for the second exchange grouping algorithm
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& projectorTree);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped second exchange communication
          */
         int mGroupedPacks2D;

      private: 
   };

   template <typename TConfigurator> BackwardSingle2DGrouper<TConfigurator>::BackwardSingle2DGrouper()
      : mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardSingle2DGrouper<TConfigurator>::~BackwardSingle2DGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardSingle2DGrouper<TConfigurator>::transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //

      // Setup the grouped communication second exchange
      this->setupGrouped2DCommunication(coord);

      //
      // Compute first step and intermediate steps of backward transform
      //
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      std::vector<Transform::TransformTree>::const_iterator it;
      for(it = coord.projectorTree().begin(); it != coord.projectorTree().end(); ++it)
      {
         // Transform scalar variable
         if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
         {
            scalIt = scalars.find(it->name());

            // Setup the first exchange communication step for scalar fields
            this->setupGrouped1DCommunication(*it, coord);

            // Compute first step of transform for scalar fields
            TConfigurator::firstStep(*it, *(scalIt->second), coord);
            // Initiate the first exchange communication step for scalar fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for scalar fields
            TConfigurator::secondStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Setup the first exchange communication step for vector fields
            this->setupGrouped1DCommunication(*it, coord);

            // Compute first step of transform for vector fields
            TConfigurator::firstStep(*it, *(vectIt->second), coord);
            // Initiate the first exchange communication step for vector fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for vector fields
            TConfigurator::secondStep(*it, *(vectIt->second), coord);
         }
      }

      // Initiate the second exchange communication step
      TConfigurator::initiate2DCommunication(coord);

      //
      // Compute last step
      //
      for(it = coord.projectorTree().begin(); it != coord.projectorTree().end(); ++it)
      {
         // Transform scalar variable
         if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
         {
            scalIt = scalars.find(it->name());

            // Compute last step of transform for scalar fields
            TConfigurator::lastStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Compute last step of transform for vector fields
            TConfigurator::lastStep(*it, *(vectIt->second), coord);
         }
      }
   }

   template <typename TConfigurator> void BackwardSingle2DGrouper<TConfigurator>::setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> id = std::make_pair(tree.name(), tree.comp<FieldComponents::Spectral::Id>());
      if(this->mNamedPacks1D.count(id) == 1)
      {
         TConfigurator::setup1DCommunication(this->mNamedPacks1D.at(id), coord);
      }
   }

   template <typename TConfigurator> void BackwardSingle2DGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks2D > 0)
      {
         TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardSingle2DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      return this->namePacks1D(projectorTree);
   }

   template <typename TConfigurator> ArrayI BackwardSingle2DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // BACKWARDSINGLE2DGROUPER_HPP

#endif //QUICC_TRANSGROUPER_SINGLE2D
