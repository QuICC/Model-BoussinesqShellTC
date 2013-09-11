/** 
 * @file BackwardTransformGrouper.hpp
 * @brief This class defines the backward transform grouping algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */
#ifdef GEOMHDISCC_TRANSGROUPER_TRANSFORM

#ifndef BACKWARDTRANSFORMGROUPER_HPP
#define BACKWARDTRANSFORMGROUPER_HPP

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
    * @brief This class defines the backward transform grouping algorithm
    */
   template <typename TConfigurator> class BackwardTransformGrouper: public IBackwardGrouper
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
          * @param varInfo Variable information
          */
         virtual ArrayI packs1D(const VariableRequirement& varInfo);

         /**
          * @brief Get the number of required packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const VariableRequirement& varInfo);

      protected:
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

      // Initiate the grouped first exchnage communication
      TConfigurator::initiate1DCommunication(coord);

      // Setup the grouped second exchange communication
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);

      //
      // Compute intermediate step
      //

      // First treat the scalar variables
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Compute first step of transform for scalar fields
         TConfigurator::secondStep(scalIt->first, *(scalIt->second), coord);
      }

      // .. then the vector variables
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
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

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI BackwardTransformGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // BACKWARDTRANSFORMGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_TRANSFORM
