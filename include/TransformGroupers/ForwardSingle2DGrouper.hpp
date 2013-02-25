/** \file ForwardSingle2DGrouper.hpp
 *  \brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
 */

#ifndef FORWARDSINGLE2DGROUPER_HPP
#define FORWARDSINGLE2DGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TransformGroupers/ForwardGrouperBase.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class ForwardSingle2DGrouper : public ForwardGrouperBase
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardSingle2DGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardSingle2DGrouper();

         /**
          * @brief Setup the full forward transform structure for the second exchange grouping algorithm
          *
          * @param scalEqs Vector of scalar equations
          * @param vectEqs Vector of vector equations
          * @param coord   Transform coord
          */
         virtual void transform(std::vector<SharedScalarEquation>& scalEqs, std::vector<SharedVectorEquation>& vectEqs, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

      protected:
         /**
          * @brief Storage for the size of the grouped second exchange communication
          */
         int mGroupedPacks2D;

      private: 
   };

   template <typename TConfigurator> ForwardSingle2DGrouper<TConfigurator>::ForwardSingle2DGrouper()
      :mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardSingle2DGrouper<TConfigurator>::~ForwardSingle2DGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardSingle2DGrouper<TConfigurator>::transform(std::vector<SharedScalarEquation>& scalEqs, std::vector<SharedVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction
      // ... and forward transform 
      //

      // Setup the grouped second exchange communication
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);

      //
      // Compute first step of forward transform
      //

      // First treat the scalar equations
      std::vector<SharedScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Compute first step of transform for scalar equation
         TConfigurator::firstStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      std::vector<SharedVectorEquation>::iterator vectEqIt;
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Compute first step of transform for vector equation
         TConfigurator::firstStep(*vectEqIt, coord);
      }

      // Initiate the grouped second exchange communication
      TConfigurator::initiate2DCommunication(coord);

      //
      // Compute intermidiate and last steps
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Setup the FDSH communication step
         TConfigurator::setup1DCommunication(this->mScalarPacks1D, coord);
         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
         // Initiate the FDSH communication step
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for scalar equation
         TConfigurator::lastStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Setup the FDSH communication step
         TConfigurator::setup1DCommunication(this->mVectorPacks1D, coord);
         // Compute second step of transform for vector equation
         TConfigurator::secondStep(*vectEqIt, coord);
         // Initiate the FDSH communication step 
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for vector equation
         TConfigurator::lastStep(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      return this->listPacks1D(varInfo);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   { 
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDSINGLE2DGROUPER_HPP
