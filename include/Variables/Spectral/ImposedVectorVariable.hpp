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

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of vector variable with an imposed component
    */
   template <typename TSScalar, typename TPScalar> class ImposedVectorVariable: public VectorVariable<TSScalar,TPScalar> 
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
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rTotal();

         /**
          * @brief Get spectral vector imposed field
          */
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  imposed() const;

         /**
          * @brief Set spectral vector imposed field 
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rImposed();

         /**
          * @brief initialise to zeros
          */
         void setZeros();

     #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
      protected:

         /**
          * @brief Spectral vector imposed field
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>    mImposed;

         /**
          * @brief Spectral vector total field
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>    mTotal;

      private:
   };

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,TPScalar>::total() const
   {
      return this->mTotal;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,TPScalar>::rTotal()
   {
      return this->mTotal;
   }

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,TPScalar>::imposed() const
   {
      return this->mImposed;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& ImposedVectorVariable<TSScalar,TPScalar>::rImposed()
   {
      return this->mImposed;
   }

   template <typename TSScalar, typename TPScalar> ImposedVectorVariable<TSScalar,TPScalar>::ImposedVectorVariable(SharedResolution spRes)
      : VectorVariable(spRes), mImposed(spRes->backwardSetup()), mTotal(spRes->backwardSetup())
   {
   }

   template <typename TSScalar, typename TPScalar> ImposedVectorVariable<TSScalar,TPScalar>::~ImposedVectorVariable()
   {
   }

   template <typename TSScalar, typename TPScalar> void ImposedVectorVariable<TSScalar,TPScalar>::setZeros()
   {
      // initialise the perturbation field to zero
      VectorVariable::setZeros();

      // initialise vector imposed field to zero
      this->rImposed().setZeros();

      // initialise vector total field to zero
      this->mTotal.setZeros();
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TSScalar, typename TPScalar> MHDFloat ImposedVectorVariable<TSScalar,TPScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += VectorVariable::requiredStorage();

      mem += this->mImposed.requiredStorage();

      mem += this->mTotal.requiredStorage();

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // IMPOSEDVECTORVARIABLE_HPP
