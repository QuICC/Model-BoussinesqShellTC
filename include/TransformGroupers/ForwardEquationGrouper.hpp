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
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      std::vector<Transform::IntegratorTree>::const_iterator it;
      for(it = coord.integratorTree().begin(); it != coord.integratorTree().end(); ++it)
      {
         // Transform scalar equation variable
         if(it->comp() == FieldComponents::Physical::SCALAR)
         {
//            scalEqIt = scalEqs.find(it->name());

            // Sychronize 
            FrameworkMacro::synchronize();

            // Setup the first exchange communication step for scalar fields
            this->setupGrouped1DCommunication(*it, coord);
            // Setup the second exchange communication step for scalar fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute first step of transform for scalar fields
            TConfigurator::firstStep(*it, *scalEqIt, coord);
            // Initiate the first exchange communication step for scalar fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for scalar fields
            TConfigurator::secondStep(*it, *scalEqIt, coord);
            // Initiate the second exchange communication step for scalar fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for scalar fields
            TConfigurator::lastStep(*it, *scalEqIt, coord);

         // Transform vector equation
         } else
         {
//            vectEqIt = vectEqs.find(it->name());

            // Sychronize 
            FrameworkMacro::synchronize();

            // Setup the first exchange communication step for vector fields
            this->setupGrouped1DCommunication(*it, coord);
            // Setup the second exchange communication step for vector fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute first step of transform for vector fields
            TConfigurator::firstStep(*it, *vectEqIt, coord);
            // Initiate the first exchange communication step for vector fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for vector fields
            TConfigurator::secondStep(*it, *vectEqIt, coord);
            // Initiate the second exchange communication step for vector fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for vector fields
            TConfigurator::lastStep(*it, *vectEqIt, coord);
         }
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

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp()));

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardEquationGrouper<TConfigurator>::setupGrouped2DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp()));

      TConfigurator::setup2DCommunication(packs, coord);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs1D(const std::vector<IntegratorTree>& integratorTree)
   {
      return this->namePacks1D(integratorTree);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs2D(const std::vector<IntegratorTree>& integratorTree)
   {  
      return this->namePacks2D(integratorTree);
   }

}
}

#endif // FORWARDEQUATIONGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_EQUATION
