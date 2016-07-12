/** 
 * @file BackwardTransformGrouper.hpp
 * @brief This class defines the backward transform grouping algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_TRANSFORM

#ifndef BACKWARDTRANSFORMGROUPER_HPP
#define BACKWARDTRANSFORMGROUPER_HPP

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

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the backward transform grouping algorithm
    */
   template <typename TConfigurator> class BackwardTransformGrouper: public IBackwardGrouper3D
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardTransformGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardTransformGrouper();

         /**
          * @brief Setup the full backward transform structure for the transform grouping algorithm
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
          * @brief Get the number of required packs for the second exchange
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
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped communication for the first exchange
          */
         int mGroupedPacks1D;

         /**
          * @brief Storage for the size of the grouped communication for the second exchange
          */
         int mGroupedPacks2D;

      private: 
   };

   template <typename TConfigurator> BackwardTransformGrouper<TConfigurator>::BackwardTransformGrouper()
      : mGroupedPacks1D(-1), mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardTransformGrouper<TConfigurator>::~BackwardTransformGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardTransformGrouper<TConfigurator>::transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //

      // Setup the 1D grouped first exchange communication
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

      // Initiate the grouped first exchange communication
      TConfigurator::initiate1DCommunication(coord);

      // Setup the grouped second exchange communication
      this->setupGrouped2DCommunication(coord);

      //
      // Compute intermediate step
      //
      for(it = coord.projectorTree().begin(); it != coord.projectorTree().end(); ++it)
      {
         // Transform scalar variable
         if(it->comp<FieldComponents::Spectral::Id>() == FieldComponents::Spectral::SCALAR)
         {
            scalIt = scalars.find(it->name());

            // Compute second step of transform for scalar fields
            TConfigurator::secondStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Compute second step of transform for vector fields
            TConfigurator::secondStep(*it, *(vectIt->second), coord);
         }
      }

      // Initiate the grouped second exchange communication
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

   template <typename TConfigurator> void BackwardTransformGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks1D > 0)
      {
         TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
      }
   }

   template <typename TConfigurator> void BackwardTransformGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      if(this->mGroupedPacks2D > 0)
      {
         TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& projectorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(projectorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // BACKWARDTRANSFORMGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_TRANSFORM
