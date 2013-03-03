/** \file BackwardSingle1DGrouper.hpp
 *  \brief This class defines the backward single grouping exchange grouping algorithm for the first exchange
 *
 *  \mhdBug Needs test
 */

#ifndef BACKWARDSINGLE1DGROUPER_HPP
#define BACKWARDSINGLE1DGROUPER_HPP

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
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);

      //
      // Compute first step of backward transform
      //

      // First treat the scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Compute first step of transform for scalar fields
         TConfigurator::firstStep(scalIt->first, *(scalIt->second), coord);
      }

      // .. then the vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
         // Compute first step of transform for vector fields
         TConfigurator::firstStep(vectIt->first, *(vectIt->second), coord);
      }

      // Initiate the grouped first exchange communication 
      TConfigurator::initiate1DCommunication(coord);

      //
      // Compute intermediate and last steps
      //

      // First treat the scalar variables
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Setup the second exchange communication step
         TConfigurator::setup2DCommunication(this->mNamedPacks2D[scalIt->first], coord);
         // Compute second step of transform for scalar fields
         TConfigurator::secondStep(scalIt->first, *(scalIt->second), coord);
         // Initiate the second exchnage communication
         TConfigurator::initiate2DCommunication(coord);

         // Compute last step of transform for scalar fields
         TConfigurator::lastStep(scalIt->first, *(scalIt->second), coord);
      }

      // .. then the vector variables
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
         // Setup the second exchange communication step
         TConfigurator::setup2DCommunication(this->mNamedPacks2D[vectIt->first], coord);
         // Compute second step of transform for vector fields
         TConfigurator::secondStep(vectIt->first, *(vectIt->second), coord);
         // Initiate the second exchange communication
         TConfigurator::initiate2DCommunication(coord);

         // Compute last step of transform for vector fields
         TConfigurator::lastStep(vectIt->first, *(vectIt->second), coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardSingle1DGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {
      return this->namePacks2D(varInfo);
   }

}
}

#endif // BACKWARDSINGLE1DGROUPER_HPP
