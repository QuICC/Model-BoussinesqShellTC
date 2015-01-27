/** 
 * @file ForwardSingle1DGrouper.hpp
 * @brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_SINGLE1D

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
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TransformGroupers/IForwardGrouper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the first exchange
    */
   template <typename TConfigurator> class ForwardSingle1DGrouper: public IForwardGrouper
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
          * @brief Setup the full forward transform structure for the first exchange grouping algorithm for the equations
          *
          * @param scalEqs Vector of scalar equations
          * @param vectEqs Vector of vector equations
          * @param coord   Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs1D(const std::vector<IntegratorTree>& integratorTree);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs2D(const std::vector<IntegratorTree>& integratorTree);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord);

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

   template <typename TConfigurator> inline void ForwardSingle1DGrouper<TConfigurator>::transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction
      // ... and forward transform 
      //
      
      // Setup the grouped first exchange communication
      this->setupGrouped1DCommunication(coord);

      //
      // Compute first and second steps of forward transform
      //

      // First treat the scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Sychronize
         FrameworkMacro::synchronize();

         bool hasNL = (*scalEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear();
         // Setup the second exchange communication step
         this->setupGrouped2DCommunication((*scalEqIt)->name(), hasNL, coord);

         // Compute first step of transform for scalar equation
         TConfigurator::firstStep(*scalEqIt, coord);
         // Initiate the second exchange communication
         TConfigurator::initiate2DCommunication(coord);

         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Sychronize
         FrameworkMacro::synchronize();

         bool hasNL = (*vectEqIt)->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear();
         // Setup the second exchange communication step
         this->setupGrouped2DCommunication((*vectEqIt)->name(), hasNL, coord);

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

      //
      // Update equation variable after transforms
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Update equation variable
         TConfigurator::updateEquation(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Update equation variable
         TConfigurator::updateEquation(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> void ForwardSingle1DGrouper<TConfigurator>::setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
   }

   template <typename TConfigurator> void ForwardSingle1DGrouper<TConfigurator>::setupGrouped2DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp()));

      TConfigurator::setup2DCommunication(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs1D(const std::vector<IntegratorTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(varInfo, nonInfo);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardSingle1DGrouper<TConfigurator>::packs2D(const std::vector<IntegratorTree>& integratorTree)
   {  
      return this->namePacks2D(integratorTree);
   }

}
}

#endif // FORWARDSINGLE1DGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_SINGLE1D
