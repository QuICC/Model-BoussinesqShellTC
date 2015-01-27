/** 
 * @file ForwardSingle2DGrouper.hpp
 * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_SINGLE2D

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
#include "TransformGroupers/IForwardGrouper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class ForwardSingle2DGrouper : public IForwardGrouper
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
          * @brief Setup the full forward transform structure for the second exchange grouping algorithm for the equations
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

   template <typename TConfigurator> inline void ForwardSingle2DGrouper<TConfigurator>::transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction
      // ... and forward transform 
      //

      // Setup the grouped second exchange communication
      this->setupGrouped2DCommunication(coord);

      //
      // Compute first step of forward transform
      //

      // First treat the scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Compute first step of transform for scalar equation
         TConfigurator::firstStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
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
         // Synchronize
         FrameworkMacro::synchronize();
         
         bool hasNL = (*scalEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear();
         // Setup the first communication step
         this->setupGrouped1DCommunication((*scalEqIt)->name(), hasNL, coord);

         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
         // Initiate the first communication step
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for scalar equation
         TConfigurator::lastStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Synchronize
         FrameworkMacro::synchronize();
         
         bool hasNL = (*vectEqIt)->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear();
         // Setup the first communication step
         this->setupGrouped1DCommunication((*vectEqIt)->name(), hasNL, coord);

         // Compute second step of transform for vector equation
         TConfigurator::secondStep(*vectEqIt, coord);
         // Initiate the first communication step 
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for vector equation
         TConfigurator::lastStep(*vectEqIt, coord);
      }

      //
      // Update equation variable after transforms
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Update equation variable after transforms for scalar equation
         TConfigurator::updateEquation(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Update equation variable after transforms for vector equation
         TConfigurator::updateEquation(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp()));

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped2DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs1D(const std::vector<IntegratorTree>& integratorTree)
   {
      return this->namePacks1D(integratorTree);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs2D(const std::vector<IntegratorTree>& integratorTree)
   { 
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo, nonInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDSINGLE2DGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_SINGLE2D
