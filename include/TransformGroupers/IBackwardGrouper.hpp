/** 
 * @file IBackwardGrouper.hpp
 * @brief This class defines some basic forward transform grouping tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBACKWARDGROUPER_HPP
#define IBACKWARDGROUPER_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/ProjectorTree.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm
    */
   class IBackwardGrouper
   {
      public:
         /// Typedef for field and component ID
         typedef std::pair<PhysicalNames::Id,FieldComponents::Spectral::Id> FieldIdType;

         /**
          * @brief Setup the full backward transform structure
          *
          * @param scalars Vector of scalar fields
          * @param vectors Vector of vector fields
          * @param coord   Transform coord
          */
         virtual void transform(std::map<PhysicalNames::Id, Datatypes::SharedScalarVariableType>& scalars, std::map<PhysicalNames::Id, Datatypes::SharedVectorVariableType>& vectors, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs1D(const std::vector<ProjectorTree>& projectorTree) = 0;

         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<ProjectorTree>& projectorTree) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Storage for named packet sizes for the first exchange
          */
         std::map<FieldIdType, int>  mNamedPacks1D;

         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<FieldIdType, int>  mNamedPacks2D;

         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks1D(const std::vector<ProjectorTree>& projectorTree);

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks2D(const std::vector<ProjectorTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks1D(const std::vector<ProjectorTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks2D(const std::vector<ProjectorTree>& projectorTree);

         /**
          * @brief Empty constructor
          */
         IBackwardGrouper();

         /**
          * @brief Empty destructor
          */
         ~IBackwardGrouper();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IBackwardGrouper>   SharedIBackwardGrouper;

}
}

#endif // IBACKWARDGROUPER_HPP
