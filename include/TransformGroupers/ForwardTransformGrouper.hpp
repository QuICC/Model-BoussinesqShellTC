/** 
 * @file ForwardTransformGrouper.hpp
 * @brief This class defines the forward transform grouping algorithm
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_TRANSFORM

#ifndef FORWARDTRANSFORMGROUPER_HPP
#define FORWARDTRANSFORMGROUPER_HPP

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
    * @brief This class defines the forward transform grouping algorithm
    */
   template <typename TConfigurator> class ForwardTransformGrouper: public IForwardGrouper
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
         void setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord);

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
      //
      // Compute nonlinear interaction
      // ... and forward transform 
      //

      // Setup the grouped second exchnage communication
      this->setupGrouped2DCommunication(coord);

      //
      // Compute first step
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

      // Setup the grouped first exchange communication
      this->setupGrouped1DCommunication(coord);

      //
      // Compute intermediate step
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Compute second step of transform for scalar equation
         TConfigurator::secondStep(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Compute second step of transform for vector equation
         TConfigurator::secondStep(*vectEqIt, coord);
      }

      // Initiate the grouped first exchange communication
      TConfigurator::initiate1DCommunication(coord);

      //
      // Compute last step
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

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped1DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);
   }

   template <typename TConfigurator> void ForwardTransformGrouper<TConfigurator>::setupGrouped2DCommunication(const IntegratorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs1D(const std::vector<IntegratorTree>& integratorTree)
   {
      // Get size of grouped communication
      ArrayI packs = this->groupPacks1D(varInfo, nonInfo);

      // Store the number of grouped packs
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs2D(const std::vector<IntegratorTree>& integratorTree)
   {  
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo, nonInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDTRANSFORMGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_TRANSFORM
