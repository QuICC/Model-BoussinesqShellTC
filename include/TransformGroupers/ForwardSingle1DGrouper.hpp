/** \file ForwardSingle1DGrouper.hpp
 *  \brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
 *
 *  \mhdBug Needs test
 */

#ifndef FORWARDSINGLE1DGROUPER_HPP
#define FORWARDSINGLE1DGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/ScalarEquation.hpp"
#include "Equations/VectorEquation.hpp"
#include "TransformGroupers/ForwardGrouperBase.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
    */
   template <typename TConfigurator> class ForwardSingle1DGrouper: public ForwardGrouperBase
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardSingle1DGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardSingle1DGrouper();

         /**
          * @brief Setup the full forward transform structure for the first exchange grouping algorithm
          *
          * @param scalEqs Vector of scalar equations
          * @param vectEqs Vector of vector equations
          * @param coord   Transform coord
          */
         virtual void transform(std::vector<SharedIScalarEquation>& scalEqs, std::vector<SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord);

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
          * @brief Storage for the size of the grouped first echange communication
          */
         int mGroupedPacks1D;

      private: 
   };

   template <typename TConfigurator> ForwardSingle1DGrouper<TConfigurator>::ForwardSingle1DGrouper()
      :mGroupedPacks1D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardSingle1DGrouper<TConfigurator>::~ForwardSingle1DGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardSingle1DGrouper<TConfigurator>::transform(std::vector<SharedIScalarEquation>& scalEqs, std::vector<SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction
      // ... and forward transform 
      //
      
      // Setup the grouped first exchange communication
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);

      //
      // Compute first and second steps of forward transform
      //

      // First treat the scalar equations
      std::vector<SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Setup the second exchange communication step
         TConfigurator::setup2DCommunication(this->mScalarPacks2D, coord);
         // Compute first step of transform for scalar equation
         TConfigurator::firstStep(*scalEqIt, coord);
         // Initiate the second exchange communication
         TConfigurator::initiate2DCommunication(coord);

         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      std::vector<SharedIVectorEquation>::iterator vectEqIt;
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Setup the second exchange communication step
         TConfigurator::setup2DCommunication(this->mVectorPacks2D, coord);
         // Compute first step of transform
         TConfigurator::firstStep(*vectEqIt, coord);
         // Initiate the second exchange communication
         TConfigurator::initiate2DCommunication(coord);

         // Compute second step of transform for vector equation
         TConfigurator::secondStep(*vectEqIt, coord);
      }

      // Initiate the grouped first exchange communication step for all equations
      TConfigurator::initiate1DCommunication(coord);

      //
      // Compute last step of forward transform
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Compute last step of transform for scalar equation
         TConfigurator::lastStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Compute last step of transform for vector equation
         TConfigurator::lastStep(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {  
      return this->listPacks2D(varInfo);
   }

}
}

#endif // FORWARDSINGLE1DGROUPER_HPP
