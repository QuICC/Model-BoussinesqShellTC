/** \file ForwardTransformGrouper.hpp
 *  \brief This class defines the forward transform grouping algorithm
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
          * @brief Setup the full forward transform structure for the transform grouping algorithm
          *
          * @param scalEqs Vector of scalar equations
          * @param vectEqs Vector of vector equations
          * @param coord   Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord);

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
      TConfigurator::setup2DCommunication(this->mGroupedPacks2D, coord);

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
      TConfigurator::setup1DCommunication(this->mGroupedPacks1D, coord);

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
      // Prepare timestep after transforms
      //

      // First treat the scalar equations
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Prepare timestep after transforms for scalar equation
         TConfigurator::prepareTimestep(*scalEqIt, coord);
      }

      // ... then the vector equations
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Prepare timestep after transforms for vector equation
         TConfigurator::prepareTimestep(*vectEqIt, coord);
      }
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      // Get size of groupe communication
      ArrayI packs = this->groupPacks1D(varInfo);

      // Store the number of grouped packs 
      this->mGroupedPacks1D = packs(0);

      return packs;
   }

   template <typename TConfigurator> ArrayI ForwardTransformGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {  
      // Get size of grouped communication
      ArrayI packs = this->groupPacks2D(varInfo);

      // Store the number of grouped packs
      this->mGroupedPacks2D = packs(0);

      return packs;
   }

}
}

#endif // FORWARDTRANSFORMGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_TRANSFORM
