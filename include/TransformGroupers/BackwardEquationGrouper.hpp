/** 
 * @file BackwardEquationGrouper.hpp
 * @brief This class defines a simple equation wise backward transform grouping algorithm (serial algorithm)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
#ifdef GEOMHDISCC_TRANSGROUPER_EQUATION

#ifndef BACKWARDDEQUATIONGROUPER_HPP
#define BACKWARDDEQUATIONGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "TransformGroupers/IBackwardGrouper.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines a simple equation wise backward transform grouping algorithm (serial algorithm)
    */
   template <typename TConfigurator> class BackwardEquationGrouper: public IBackwardGrouper
   {
      public:
         /**
          * @brief Constructor
          */
         BackwardEquationGrouper();

         /**
          * @brief Destructor
          */
         ~BackwardEquationGrouper();

         /**
          * @brief Setup the full backward transform structure
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord);

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

      private: 
   };

   template <typename TConfigurator> BackwardEquationGrouper<TConfigurator>::BackwardEquationGrouper()
   {
      this->split = TConfigurator::SplitLocation;
   }

   template <typename TConfigurator> BackwardEquationGrouper<TConfigurator>::~BackwardEquationGrouper()
   {
   }

   template <typename TConfigurator> inline void BackwardEquationGrouper<TConfigurator>::transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord)
   {
      //
      // Compute backward transform
      //
      
      // First treat the scalar variables
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      for(scalIt = scalars.begin(); scalIt != scalars.end(); scalIt++)
      {
         // Setup the first exchange communication step for scalar fields
         TConfigurator::setup1DCommunication(this->mNamedPacks1D.at(scalIt->first), coord);
         // Setup the second exchange communication step for scalar fields
         TConfigurator::setup2DCommunication(this->mNamedPacks2D.at(scalIt->first), coord);

         // Compute first step of transform for scalar fields
         TConfigurator::firstStep(scalIt->first, *(scalIt->second), coord);
         // Initiate the first exchange communication step for scalar fields
         TConfigurator::initiate1DCommunication(coord);

         // Compute second step of transform for scalar fields
         TConfigurator::secondStep(scalIt->first, *(scalIt->second), coord);
         // Initiate the second exchange communication step for scalar fields
         TConfigurator::initiate2DCommunication(coord);

         // Compute last step of transform for scalar fields
         TConfigurator::lastStep(scalIt->first, *(scalIt->second), coord);
      }
      
      // .. then the vector variables
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      for(vectIt = vectors.begin(); vectIt != vectors.end(); vectIt++)
      {
         // Setup the first exchange communication step for vector fields
         TConfigurator::setup1DCommunication(this->mNamedPacks1D.at(vectIt->first), coord);
         // Setup the second exchange communication step for vector fields
         TConfigurator::setup2DCommunication(this->mNamedPacks2D.at(vectIt->first), coord);

         // Compute first step of transform for vector fields
         TConfigurator::firstStep(vectIt->first, *(vectIt->second), coord);
         // Initiate the first exchange communication step for vector fields
         TConfigurator::initiate1DCommunication(coord);

         // Compute second step of transform for vector fields
         TConfigurator::secondStep(vectIt->first, *(vectIt->second), coord);
         // Initiate the second exchange communication step for vector fields
         TConfigurator::initiate2DCommunication(coord);

         // Compute last step of transform for vector fields
         TConfigurator::lastStep(vectIt->first, *(vectIt->second), coord);
      }
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs1D(const VariableRequirement& varInfo)
   {
      return this->namePacks1D(varInfo);
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs2D(const VariableRequirement& varInfo)
   {
      return this->namePacks2D(varInfo);
   }

}
}

#endif // BACKWARDDEQUATIONGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_EQUATION
