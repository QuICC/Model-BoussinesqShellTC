/** 
 * @file ImposedVectorVariable.hpp
 * @brief Implementation of vector variable with an imposed component
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IMPOSEDVECTORVARIABLE_HPP
#define IMPOSEDVECTORVARIABLE_HPP

// Configuration includes
//
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"
#include "Variables/Spectral/VectorVariable.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Implementation of vector variable with an imposed component
    */
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> class ImposedVectorVariable: public VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS> 
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         ImposedVectorVariable(SharedReslolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~ImposedVectorVariable();

         /**
          * @brief Get spectral vector field (total field)
          */
         const VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>&  rTotal();

         /**
          * @brief Get spectral vector imposed field
          */
         const VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>&  imposed() const;

         /**
          * @brief Set spectral vector imposed field 
          */
         VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>&  rImposed();

         /**
          * @brief initialise to zeros
          */
         void setZeros();

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

         /**
          * @brief Spectral vector imposed field
          */
         VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>    mImposed;

         /**
          * @brief Spectral vector total field
          */
         VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>    mTotal;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::total() const
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rTotal()
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::imposed() const
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rImposed()
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::ImposedVectorVariable(SharedResolution spRes)
      : VectorVariable(spRes), mImposed(spRes->backwardSetup()), mTotal(spRes->backwardSetup())
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::~ImposedVectorVariable()
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::setZeros()
   {
      // initialise the perturbation field to zero
      VectorVariable::setZeros();

      // initialise vector imposed field to zero
      this->rImposed().setZeros();

      // initialise vector total field to zero
      this->mTotal.setZeros();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> MHDFloat ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += VectorVariable::requiredStorage();

      mem += this->mImposed.requiredStorage();

      mem += this->mTotal.requiredStorage();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // IMPOSEDVECTORVARIABLE_HPP
