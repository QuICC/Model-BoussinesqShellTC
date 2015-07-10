/** 
 * @file IForwardGrouper2D.hpp
 * @brief This class defines some basic forward transform grouping tools in 2D Space
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IFORWARDGROUPER2D_HPP
#define IFORWARDGROUPER2D_HPP

// Configuration includes
// 
#include "TypeSelectors/TransformCommSelector.hpp"
#include "TypeSelectors/TransformTreeSelector.hpp"

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/IScalarEquation.hpp"
#include "Equations/IVectorEquation.hpp"

namespace GeoMHDiSCC {

namespace Transform {

   /**
    * @brief This class defines the serial forward transform grouping algorithm in 2D space
    */
   class IForwardGrouper2D
   {
      public:
         /// Typedef for field and component ID
         typedef std::pair<PhysicalNames::Id,FieldComponents::Physical::Id> FieldIdType;

         /**
          * @brief Setup the full forward transform structure for equations
          *
          * @param scalEqs  Vector of scalar equations
          * @param vectEqs  Vector of vector equations
          * @param coord    Transform coord
          */
         virtual void transform(std::vector<Equations::SharedIScalarEquation>& scalEqs, std::vector<Equations::SharedIVectorEquation>& vectEqs, TransformCoordinatorType& coord) = 0;

         /**
          * @brief Get the number of required buffer packs for the first exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         virtual ArrayI packs1D(const std::vector<IntegratorTree>& integratorTree) = 0;

         /**
          * @brief Location of the split in the configurator
          */
         Splitting::Locations::Id split;

      protected:
         /**
          * @brief Find equation corresponding to name of scalar
          */
         template <typename TSharedEquation> typename std::vector<TSharedEquation>::iterator findEquation(std::vector<TSharedEquation>& eqs, PhysicalNames::Id name);

         /**
          * @brief Storage for named packet sizes for the first exchange
          */
         std::map<FieldIdType, int>  mNamedPacks1D;

         /**
          * @brief Get and set the name pack numbers for the first exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         ArrayI namePacks1D(const std::vector<IntegratorTree>& integratorTree);

         /**
          * @brief Get the grouped pack number for the first exchange
          *
          * @param integratorTree Transform integrator tree 
          */
         ArrayI groupPacks1D(const std::vector<IntegratorTree>& integratorTree);

         /**
          * @brief Empty constructor
          */
         IForwardGrouper2D();

         /**
          * @brief Empty destructor
          */
         ~IForwardGrouper2D();

      private: 
   };

   template < typename TSharedEquation> typename std::vector<TSharedEquation>::iterator IForwardGrouper2D::findEquation(std::vector<TSharedEquation>& eqs, PhysicalNames::Id name)
   {
      typename std::vector<TSharedEquation>::iterator eqIt;
      for(eqIt = eqs.begin(); eqIt != eqs.end(); ++eqIt)
      {
         if((*eqIt)->name() == name)
         {
            break;
         }
      }

      return eqIt;
   }

   /// Typdef for a smart reference counting pointer to a backward grouper base
   typedef SharedPtrMacro<IForwardGrouper2D>   SharedIForwardGrouper2D;

   #ifdef GEOMHDISCC_SPATIALDIMENSION_2D
      typedef IForwardGrouper2D IForwardGrouper;

      typedef SharedPtrMacro<IForwardGrouper2D>   SharedIForwardGrouper;
   #endif //GEOMHDISCC_SPATIALDIMENSION_2D

}
}

#endif // IFORWARDGROUPER2D_HPP
