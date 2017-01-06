/** 
 * @file BackwardSingle1DGrouper.hpp
 * @brief This class defines the backward single grouping exchange grouping algorithm for the first exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef QUICC_TRANSGROUPER_SINGLE1D

#ifndef BACKWARDSINGLE1DGROUPER_HPP
#define BACKWARDSINGLE1DGROUPER_HPP

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
#include "TransformGroupers/IBackwardGrouperMacro.h"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the backward single grouping exchange grouping algorithm for the first exchange
    */
   template <typename TConfigurator> class BackwardSingle1DGrouper: public IBackwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardSingle1DGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardSingle1DGrouper();

         /**
          * @brief Setup the full backward transform structure for the first exchange grouping algorithm
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
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped first exchange communication
          */
         int mGroupedPacks1D;

      private: 
   };

   template <typename TConfigurator> BackwardSingle1DGrouper<TConfigurator>::BackwardSingle1DGrouper()
      : mGroupedPacks1D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardSingle1DGrouper<TConfigurator>::~BackwardSingle1DGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardSingle1DGrouper<TConfigurator>::transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //

      // Setup the grouped first exchange communication
      this->setupGrouped1DCommunication(coord);

      //
      // Compute first step of backward transform
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

            // Compute first step of transform for scalar fields
            TConfigurator::firstStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Compute first step of transform for vector fields
            TConfigurator::firstStep(*it, *(vectIt->second), coord);
         }
      }

      // Initiate the first exchange communication step for vector fields
      TConfigurator::initiate1DCommunication(coord);

      //
      // Compute intermediate and last steps
      //
      for(it = coord.projectorTree().begin(); it != coord.projectorTree().end(); ++it)
      {
         // Transform scalar variable
         if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
         {
            scalIt = scalars.find(it->name());

            // Setup the second exchange communication step for scalar fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute second step of transform for scalar fields
            TConfigurator::secondStep(*it, *(scalIt->second), coord);
            // Initiate the second exchange communication step for scalar fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for scalar fields
            TConfigurator::lastStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Setup the second exchange communication step for vector fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute second step of transform for vector fields
            TConfigurator::secondStep(*it, *(vectIt->second), coord);
            // Initiate the second exchange communication step for vector fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for vector fields
            TConfigurator::lastStep(*it, *(vectIt->second), coord);
         }
      }
   }

   template <typename TConfigurator> void BackwardSingle1DGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks1D > 0)
      {
         TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
      }
   }

   template <typename TConfigurator> void BackwardSingle1DGrouper<TConfigurator>::setupGrouped2DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> id = std::make_pair(tree.name(), tree.comp<FieldComponents::Spectral::Id>());
      if(this->mNamedPacks2D.count(id) == 1)
      {
         TConfigurator::setup2DCommunication(this->mNamedPacks2D.at(id), coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      return this->namePacks2D(projectorTree);
   }

}
}

#endif // BACKWARDSINGLE1DGROUPER_HPP

#endif //QUICC_TRANSGROUPER_SINGLE1D
