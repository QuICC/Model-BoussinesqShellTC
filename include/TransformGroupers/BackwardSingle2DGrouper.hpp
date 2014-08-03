/** 
 * @file BackwardSingle2DGrouper.hpp
 * @brief This class defines the backward single grouping exchange grouping algorithm for the second exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_SINGLE2D

#ifndef BACKWARDSINGLE2DGROUPER_HPP
#define BACKWARDSINGLE2DGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "TransformGroupers/IBackwardGrouper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the backward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class BackwardSingle2DGrouper: public IBackwardGrouper
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
          * @param varInfo Variable information
          */
         virtual ArrayI packs1D(const VariableRequirement& varInfo);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const VariableRequirement& varInfo);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const PhysicalNames::Id id, TransformCoordinatorType& coord);

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

      // First treat the scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Synchronize
         FrameworkMacro::synchronize();
         
         // Setup the first exchange communication
         this->setupGrouped1DCommunication(scalIt->first, coord);
         // Compute first step of transform for scalar fields
         TConfigurator::firstStep(scalIt->first, *(scalIt->second), coord);
         // Initiate the first exchange communication
         TConfigurator::initiate1DCommunication(coord);

         // Compute second step of transform for scalar fields
         TConfigurator::secondStep(scalIt->first, *(scalIt->second), coord);
      }

      // .. then the vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
         // Synchronize
         FrameworkMacro::synchronize();
         
         // Setup the first exchange communication
         this->setupGrouped1DCommunication(vectIt->first, coord);
         // Compute first step of transform for vector fields
         TConfigurator::firstStep(vectIt->first, *(vectIt->second), coord);
         // Initiate the first exchange communication
         TConfigurator::initiate1DCommunication(coord);

         // Compute second step of transform for vector fields
         TConfigurator::secondStep(vectIt->first, *(vectIt->second), coord);
      }

      // Initiate the grouped second exchange communication
      TConfigurator::initiate2DCommunication(coord);

      //
      // Compute last step
      //

      // First treat the scalar variables
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Compute last step of transform for scalar fields
         TConfigurator::lastStep(scalIt->first, *(scalIt->second), coord);
      }

      // .. then the vector variables
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
         // Compute last step of transform for vector fields
         TConfigurator::lastStep(vectIt->first, *(vectIt->second), coord);
      }
   }

   template <typename TConfigurator> void BackwardSingle2DGrouper<TConfigurator>::setupGrouped1DCommunication(const PhysicalNames::Id id, TransformCoordinatorType& coord)
   {
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

   template <typename TConfigurator> ArrayI BackwardSingle2DGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      return this->namePacks1D(varInfo);
   }

   template <typename TConfigurator> ArrayI BackwardSingle2DGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // BACKWARDSINGLE2DGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_SINGLE2D
