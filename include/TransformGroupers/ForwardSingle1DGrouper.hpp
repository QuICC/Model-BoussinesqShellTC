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
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

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
      // Setup the second exchange communication step for scalar fields
      this->setupGrouped1DCommunication(coord);

      //
      // Compute nonlinear interaction 
      // ... and firts and second forward transform steps
      //
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      std::vector<Transform::IntegratorTree>::const_iterator it;
      for(it = coord.integratorTree().begin(); it != coord.integratorTree().end(); ++it)
      {
         // Transform scalar equation variable
         if(it->comp() == FieldComponents::Physical::SCALAR)
         {
            scalEqIt = this->findEquation(scalEqs,it->name());

            // Presolve might not provide all equations
            if(scalEqIt != scalEqs.end())
            {
               // Sychronize 
               FrameworkMacro::synchronize();

               // Setup the first exchange communication step for scalar fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, *scalEqIt, coord);
               // Initiate the first exchange communication step for scalar fields
               TConfigurator::initiate2DCommunication(coord);

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, *scalEqIt, coord);
            }

         // Transform vector equation
         } else
         {
            vectEqIt = this->findEquation(vectEqs,it->name());

            // Presolve might not provide all equations
            if(vectEqIt != vectEqs.end())
            {
               // Sychronize 
               FrameworkMacro::synchronize();

               // Setup the first exchange communication step for vector fields
               this->setupGrouped2DCommunication(*it, coord);

               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, *vectEqIt, coord);
               // Initiate the first exchange communication step for vector fields
               TConfigurator::initiate2DCommunication(coord);

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, *vectEqIt, coord);
            }
         }
      }

      // Initiate the second exchange communication step for scalar fields
      TConfigurator::initiate1DCommunication(coord);

      //
      // ... and last step of forward transform
      //
      for(it = coord.integratorTree().begin(); it != coord.integratorTree().end(); ++it)
      {
         // Transform scalar equation variable
         if(it->comp() == FieldComponents::Physical::SCALAR)
         {
            scalEqIt = this->findEquation(scalEqs,it->name());

            // Presolve might not provide all equations
            if(scalEqIt != scalEqs.end())
            {
               // Compute last step of transform for scalar fields
               TConfigurator::lastStep(*it, *scalEqIt, coord);
            }

         // Transform vector equation
         } else
         {
            vectEqIt = this->findEquation(vectEqs,it->name());

            // Presolve might not provide all equations
            if(vectEqIt != vectEqs.end())
            {
               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, *vectEqIt, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void ForwardSingle1DGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
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
      ArrayI packs = this->groupPacks1D(integratorTree);

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
