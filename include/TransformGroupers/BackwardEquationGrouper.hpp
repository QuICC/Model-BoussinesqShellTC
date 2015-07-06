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
#include "Framework/FrameworkMacro.h"
#include "TypeSelectors/TransformCommSelector.hpp"
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
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs1D(const std::vector<ProjectorTree>& projectorTree);

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<ProjectorTree>& projectorTree);

      protected:
         /**
          * @brief Setup grouped first exchange communication
          */
         void setupGrouped1DCommunication(const ProjectorTree& tree, TransformCoordinatorType& coord);

         /**
          * @brief Setup grouped second exchange communication
          */
         void setupGrouped2DCommunication(const ProjectorTree& tree, TransformCoordinatorType& coord);

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
      std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>::iterator scalIt;
      std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>::iterator vectIt;
      std::vector<Transform::ProjectorTree>::const_iterator it;
      for(it = coord.projectorTree().begin(); it != coord.projectorTree().end(); ++it)
      {
         // Transform scalar variable
         if(it->comp() == FieldComponents::Spectral::SCALAR)
         {
            scalIt = scalars.find(it->name());

            // Sychronize 
            FrameworkMacro::synchronize();

            // Setup the first exchange communication step for scalar fields
            this->setupGrouped1DCommunication(*it, coord);
            // Setup the second exchange communication step for scalar fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute first step of transform for scalar fields
            TConfigurator::firstStep(*it, *(scalIt->second), coord);
            // Initiate the first exchange communication step for scalar fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for scalar fields
            TConfigurator::secondStep(*it, *(scalIt->second), coord);
            // Initiate the second exchange communication step for scalar fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for scalar fields
            TConfigurator::lastStep(*it, *(scalIt->second), coord);

         // Transform vector variable
         } else
         {
            vectIt = vectors.find(it->name());

            // Sychronize 
            FrameworkMacro::synchronize();

            // Setup the first exchange communication step for vector fields
            this->setupGrouped1DCommunication(*it, coord);
            // Setup the second exchange communication step for vector fields
            this->setupGrouped2DCommunication(*it, coord);

            // Compute first step of transform for vector fields
            TConfigurator::firstStep(*it, *(vectIt->second), coord);
            // Initiate the first exchange communication step for vector fields
            TConfigurator::initiate1DCommunication(coord);

            // Compute second step of transform for vector fields
            TConfigurator::secondStep(*it, *(vectIt->second), coord);
            // Initiate the second exchange communication step for vector fields
            TConfigurator::initiate2DCommunication(coord);

            // Compute last step of transform for vector fields
            TConfigurator::lastStep(*it, *(vectIt->second), coord);
         }
      }
   }

   template <typename TConfigurator> void BackwardEquationGrouper<TConfigurator>::setupGrouped1DCommunication(const ProjectorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup1DCommunication(this->mNamedPacks1D.at(std::make_pair(tree.name(), tree.comp())), coord);
   }

   template <typename TConfigurator> void BackwardEquationGrouper<TConfigurator>::setupGrouped2DCommunication(const ProjectorTree& tree, TransformCoordinatorType& coord)
   {
      TConfigurator::setup2DCommunication(this->mNamedPacks2D.at(std::make_pair(tree.name(), tree.comp())), coord);
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs1D(const std::vector<ProjectorTree>& projectorTree)
   {
      return this->namePacks1D(projectorTree);
   }

   template <typename TConfigurator> ArrayI BackwardEquationGrouper<TConfigurator>::packs2D(const std::vector<ProjectorTree>& projectorTree)
   {
      return this->namePacks2D(projectorTree);
   }

}
}

#endif // BACKWARDDEQUATIONGROUPER_HPP

#endif //GEOMHDISCC_TRANSGROUPER_EQUATION
