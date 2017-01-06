/** 
 * @file ForwardTransformGrouper.hpp
 * @brief This class defines the forward transform grouping algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef QUICC_TRANSGROUPER_TRANSFORM

#ifndef FORWARDTRANSFORMGROUPER_HPP
#define FORWARDTRANSFORMGROUPER_HPP

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

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the forward transform grouping algorithm
    */
   template <typename TConfigurator> class ForwardTransformGrouper: public IForwardGrouper3D
   {
      public:
         /**
          * @brief Constructor
          */
         ForwardTransformGrouper();

         /**
          * @brief Destructor
          */
         ~ForwardTransformGrouper();

         /**
          * @brief Setup the full forward transform structure for the transform grouping algorithm for the equations
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
         void setupGrouped1DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(TransformCoordinatorType& coord);

         /**
          * @brief Storage for the size of the grouped first exchange communication
          */
         int mGroupedPacks1D;

         /**
          * @brief Storage for the size of the grouped second exchange communication
          */
         int mGroupedPacks2D;

      private: 
   };

   template <typename TConfigurator> ForwardTransformGrouper<TConfigurator>::ForwardTransformGrouper()
      : mGroupedPacks1D(-1), mGroupedPacks2D(-1)
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardTransformGrouper<TConfigurator>::~ForwardTransformGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardTransformGrouper<TConfigurator>::transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord)
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

      // Initiate the first exchange communication step for scalar fields
      TConfigurator::initiate2DCommunication(coord);

      // Setup the second exchange communication step for scalar fields
      this->setupGrouped1DCommunication(coord);

      //
      // ... and second step of forward transform
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
               // Compute second step of transform for vector fields
               TConfigurator::secondStep(*it, *vectEqIt, coord);
            }
         }
      }

      // Initiate the second exchange communication step for vector fields
      TConfigurator::initiate1DCommunication(coord);

      //
      // ... and last step of forward transform
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

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped1DCommunication(TransformCoordinatorType& coord)
   {
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
   }

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped2DCommunication(TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs1D(const std::vector<TransformTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs2D(const std::vector<TransformTree>& integratorTree)
   {  
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(integratorTree);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDTRANSFORMGROUPER_HPP

#endif //QUICC_TRANSGROUPER_TRANSFORM
