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
#include "TypeSelectors/TransformCommSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"
#include "TransformGroupers/IForwardGrouper3D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the forward single grouping exchange grouping algorithm for the second exchange
    */
   template <typename TConfigurator> class ForwardSingle2DGrouper : public IForwardGrouper3D
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
         virtual ArrayI packs1D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& integratorTree);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord);

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
      // Setup the first exchange communication step for scalar fields
      this->setupGrouped2DCommunication(coord);

      //
      // Compute nonlinear interaction 
      // ... and first step of forward transform
      //
      std::vector<Equations::SharedIScalarEquation>::iterator scalEqIt;
      std::vector<Equations::SharedIVectorEquation>::iterator vectEqIt;
      std::vector<Transform::TransformTree>::const_iterator it;
      for(it = coord.integratorTree().begin(); it != coord.integratorTree().end(); ++it)
      {
         // Transform scalar equation variable
         if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
         {
            scalEqIt = this->findEquation(scalEqs,it->name());

            // Presolve might not provide all equations
            if(scalEqIt != scalEqs.end())
            {
               // Compute first step of transform for scalar fields
               TConfigurator::firstStep(*it, *scalEqIt, coord);
            }

         // Transform vector equation
         } else
         {
            vectEqIt = this->findEquation(vectEqs,it->name());

            // Presolve might not provide all equations
            if(vectEqIt != vectEqs.end())
            {
               // Compute first step of transform for vector fields
               TConfigurator::firstStep(*it, *vectEqIt, coord);
            }
         }
      }

      // Initiate the first exchange communication step for vector fields
      TConfigurator::initiate2DCommunication(coord);

      //
      // ... and second and last step of forward transform
      //
      for(it = coord.integratorTree().begin(); it != coord.integratorTree().end(); ++it)
      {
         // Transform scalar equation variable
         if(it->comp<FieldComponents::Physical::Id>() == FieldComponents::Physical::SCALAR)
         {
            scalEqIt = this->findEquation(scalEqs,it->name());

            // Presolve might not provide all equations
            if(scalEqIt != scalEqs.end())
            {
               // Setup the second exchange communication step for scalar fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute second step of transform for scalar fields
               TConfigurator::secondStep(*it, *scalEqIt, coord);
               // Initiate the second exchange communication step for scalar fields
               TConfigurator::initiate1DCommunication(coord);

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
               // Setup the second exchange communication step for vector fields
               this->setupGrouped1DCommunication(*it, coord);

               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, *vectEqIt, coord);
               // Initiate the second exchange communication step for vector fields
               TConfigurator::initiate1DCommunication(coord);

               // Compute last step of transform for vector fields
               TConfigurator::lastStep(*it, *vectEqIt, coord);
            }
         }
      }
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped1DCommunication(const TransformTree& tree, TransformCoordinatorType& coord)
   {
      int packs = this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp<FieldComponents::Physical::Id>()));

      TConfigurator::setup1DCommunication(packs, coord);
   }

   template <typename TConfigurator> void ForwardSingle2DGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      return this->namePacks1D(integratorTree);
   }

   template <typename TConfigurator> ArrayI ForwardSingle2DGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   { 
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDSINGLE2DGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_SINGLE2D
