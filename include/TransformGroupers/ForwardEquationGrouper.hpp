/** 
 * @file ForwardEquationGrouper.hpp
 * @brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_EQUATION

#ifndef FORWARDEQUATIONGROUPER_HPP
#define FORWARDEQUATIONGROUPER_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
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
    * @brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
    */
   template <typename TConfigurator> class ForwardEquationGrouper: public IForwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardEquationGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardEquationGrouper();

         /**
          * @brief Setup the full forward transform structure for the equations
          *
          * @param scalEqs  Vector of scalar equations
          * @param vectEqs  Vector of vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         virtual ArrayI packs1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          * @param nonInfo Nonlinear requirements
          */
         virtual ArrayI packs2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const PhysicalNames::Id id, const bool hasNL, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const PhysicalNames::Id id, const bool hasNL, TransformCoordinatorType& coord);

      private: 
   };

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::ForwardEquationGrouper()
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::~ForwardEquationGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardEquationGrouper<TConfigurator>::transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction 
      // ... and forward transform
      //

      // First treat the scalar equations
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt != scalEqs.end(); scalEqIt++)
      {  
         // Sychronize 
         FrameworkMacro::synchronize();

         bool hasNL = (*scalEqIt)->couplingInfo(FieldComponents::Spectral::SCALAR).hasNonlinear();
        // Setup the second exchange communication step for scalar equation
         this->setupGrouped2DCommunication((*scalEqIt)->name(), hasNL, coord);
         // Setup the first exchange communication step for scalar equation
         this->setupGrouped1DCommunication((*scalEqIt)->name(), hasNL, coord);

         // Compute first step of transform for scalar equation
         TConfigurator::firstStep(*scalEqIt, coord);
         // Initiate the second exchange communication step for scalar equation
         TConfigurator::initiate2DCommunication(coord);

         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
         // Initiate the first exchange communication step for scalar equation
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for scalar equation
         TConfigurator::lastStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      for(vectEqIt = vectEqs.begin(); vectEqIt != vectEqs.end(); vectEqIt++)
      {
         // Sychronize 
         FrameworkMacro::synchronize();

         bool hasNL = (*vectEqIt)->couplingInfo(FieldComponents::Spectral::ONE).hasNonlinear();
         // Setup the second exchange communication step for vector equation
         this->setupGrouped2DCommunication((*vectEqIt)->name(), hasNL, coord);
         // Setup the first exchange communication step for vector equation
         this->setupGrouped1DCommunication((*vectEqIt)->name(), hasNL, coord);

         // Compute first step of transform for vector equation
         TConfigurator::firstStep(*vectEqIt, coord);
         // Initiate the second exchange communication step for the vector equation
         TConfigurator::initiate2DCommunication(coord);

         // Compute second step of transform for vector equation
         TConfigurator::secondStep(*vectEqIt, coord);
         // Initiate the first exchange communication step for the vector equation
         TConfigurator::initiate1DCommunication(coord);

         // Compute last step of transform for vector equation
         TConfigurator::lastStep(*vectEqIt, coord);
      }

      //
      // Update equation variable after transforms
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt != scalEqs.end(); scalEqIt++)
      {
         // Update equation variable after transforms for scalar equation
         TConfigurator::updateEquation(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt != vectEqs.end(); vectEqIt++)
      {
         // Update equation variable after transforms for vector equation
         TConfigurator::updateEquation(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped1DCommunication(const PhysicalNames::Id id, const bool hasNL, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(id);
      if(!hasNL)
      {
         packs = 0;
      }

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped2DCommunication(const PhysicalNames::Id id, const bool hasNL, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks2D.at(id);
      if(!hasNL)
      {
         packs = 0;
      }

      TConfigurator::setup2DCommunication(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {
      return this->namePacks1D(varInfo, nonInfo);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo, const std::set<PhysicalNames::Id>& nonInfo)
   {  
      return this->namePacks2D(varInfo, nonInfo);
   }

}
}

#endif // FORWARDEQUATIONGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_EQUATION
