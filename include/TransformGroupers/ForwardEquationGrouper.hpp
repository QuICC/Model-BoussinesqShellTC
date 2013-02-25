/** \file ForwardEquationGrouper.hpp
 *  \brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
 */

#ifndef FORWARDEQUATIONGROUPER_HPP
#define FORWARDEQUATIONGROUPER_HPP

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
    * @brief This class defines a simple equation wise forward transform grouping algorithm (serial algorithm)
    */
   template <typename TConfigurator> class ForwardEquationGrouper: public ForwardGrouperBase
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
          * @brief Setup the full forward transform structure
          *
          * @param scalEqs  Vector of scalar equations
          * @param vectEqs  Vector of vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<SharedScalarEquation>& scalEqs, std::vector<SharedVectorEquation>& vectEqs, TransformCoordinatorType& coord);

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param varInfo Variable information
          */
         virtual ArrayI packs2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo);

      protected:

      private: 
   };

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::ForwardEquationGrouper()
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> ForwardEquationGrouper<TConfigurator>::~ForwardEquationGrouper()
   {
   }

   template <typename TConfigurator> inline void ForwardEquationGrouper<TConfigurator>::transform(std::vector<SharedScalarEquation>& scalEqs, std::vector<SharedVectorEquation>& vectEqs, TransformCoordinatorType& coord)
   {
      //
      // Compute nonlinear interaction 
      // ... and forward transform
      // ... and preparation of timestep
      //

      // First treat the scalar equations
      std::vector<SharedScalarEquation>::iterator scalEqIt;
      for(scalEqIt = scalEqs.begin(); scalEqIt < scalEqs.end(); scalEqIt++)
      {
         // Setup the second exchange communication step for scalar equation
         TConfigurator::setup2DCommunication(this->mScalarPacks2D, coord);
         // Setup the first exchange communication step for scalar equation
         TConfigurator::setup1DCommunication(this->mScalarPacks1D, coord);

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
      std::vector<SharedVectorEquation>::iterator vectEqIt;
      for(vectEqIt = vectEqs.begin(); vectEqIt < vectEqs.end(); vectEqIt++)
      {
         // Setup the second exchange communication step for vector equation
         TConfigurator::setup2DCommunication(this->mVectorPacks2D, coord);
         // Setup the first exchange communication step for vector equation
         TConfigurator::setup1DCommunication(this->mVectorPacks1D, coord);

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
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs1D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {
      return this->listPacks1D(varInfo);
   }

   template <typename TConfigurator> ArrayI ForwardEquationGrouper<TConfigurator>::packs2D(const std::map<PhysicalNames::Id, std::pair<bool,TriBool> >& varInfo)
   {  
      return this->listPacks2D(varInfo);
   }

}
}

#endif // FORWARDEQUATIONGROUPER_HPP
