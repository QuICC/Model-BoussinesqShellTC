/** 
 * @file IForwardGrouper3D.hpp
 * @brief This class defines some basic forward transform grouping tools in 3D Space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IFORWARDGROUPER3D_HPP
#define IFORWARDGROUPER3D_HPP

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
#include "TransformConfigurators/TransformTree.hpp"
#include "TransformGroupers/IForwardGrouper2D.hpp"

namespace QuICC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in 3D space
    */
   class IForwardGrouper3D: public IForwardGrouper2D
   {
      public:
         /**
          * @brief Get the number of required buffer packs for the second exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         virtual ArrayI packs2D(const std::vector<TransformTree>& integratorTree) = 0;

      protected:
         /**
          * @brief Storage for named packet sizes for the second exchange
          */
         std::map<FieldIdType, int>  mNamedPacks2D;

         /**
          * @brief Get and set the named pack numbers for the second exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         ArrayI namePacks2D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Get the grouped pack number for the second exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         ArrayI groupPacks2D(const std::vector<TransformTree>& integratorTree);

         /**
          * @brief Empty constructor
          */
         IForwardGrouper3D();

         /**
          * @brief Empty destructor
          */
         ~IForwardGrouper3D();

      private: 
   };

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IForwardGrouper3D>   SharedIForwardGrouper3D;

   #ifdef QUICC_SPATIALDIMENSION_3D
      typedef IForwardGrouper3D IForwardGrouper;

      typedef SharedPtrMacro<IForwardGrouper3D>   SharedIForwardGrouper;
   #endif //QUICC_SPATIALDIMENSION_3D

}
}

#endif // IFORWARDGROUPER3D_HPP
