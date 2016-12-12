/** 
 * @file IBackwardGrouper3D.hpp
 * @brief This class defines some basic forward transform grouping tools in 3D space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBACKWARDGROUPER3D_HPP
#define IBACKWARDGROUPER3D_HPP

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
#include "TransformGroupers/IBackwardGrouper2D.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in 3D space
    */
   class IBackwardGrouper3D: public IBackwardGrouper2D
   {
      public:
         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& projectorTree) = 0;

      protected:
         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<FieldIdType, int>  mNamedPacks2D;

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI namePacks2D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param projectorTree Transform projector tree
          */
         ArrayI groupPacks2D(const std::vector<TransformTree>& projectorTree);

         /**
          * @brief Empty constructor
          */
         IBackwardGrouper3D();

         /**
          * @brief Empty destructor
          */
         ~IBackwardGrouper3D();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IBackwardGrouper3D>   SharedIBackwardGrouper3D;

   #ifdef GEOMHDISCC_SPATIALDIMENSION_3D
      typedef IBackwardGrouper3D IBackwardGrouper;

      typedef SharedPtrMacro<IBackwardGrouper3D>   SharedIBackwardGrouper;
   #endif //GEOMHDISCC_SPATIALDIMENSION_3D

}
}

#endif // IBACKWARDGROUPER3D_HPP
