/** 
 * @file IBackwardGrouper2D.hpp
 * @brief This class defines some basic forward transform grouping tools in 2D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBACKWARDGROUPER2D_HPP
#define IBACKWARDGROUPER2D_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "TransformConfigurators/TransformTree.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in 2D space
    */
   class IBackwardGrouper2D
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
         virtual ArrayI packs1D(const std::vector<TransformTree>& projectorTree) = 0;

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
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks1D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Empty constructor
          */
         IBackwardGrouper2D();

         /**
          * @brief Empty destructor
          */
         ~IBackwardGrouper2D();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IBackwardGrouper2D>   SharedIBackwardGrouper2D;

   #ifdef QUICC_SPATIALDIMENSION_2D
      typedef IBackwardGrouper2D IBackwardGrouper;

      typedef SharedPtrMacro<IBackwardGrouper2D>   SharedIBackwardGrouper;
   #endif //QUICC_SPATIALDIMENSION_2D

}
}

#endif // IBACKWARDGROUPER2D_HPP
