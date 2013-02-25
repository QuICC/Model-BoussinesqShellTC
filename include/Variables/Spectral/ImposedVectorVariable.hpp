/** \file ImposedVectorVariable.hpp
 *  \brief Implementation of vector variable with an imposed component
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

namespace Dataypes {

   /**
    * \brief Implementation of vector variable with an imposed component
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
         const VectorField<TSScalar,SCOMPONENTS>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,SCOMPONENTS>&  rTotal();

         /**
          * @brief Get spectral vector imposed field
          */
         const VectorField<TSScalar,SCOMPONENTS>&  imposed() const;

         /**
          * @brief Set spectral vector imposed field 
          */
         VectorField<TSScalar,SCOMPONENTS>&  rImposed();

         /**
          * @brief initialise to zeros
          */
         void initialiseZeros();

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
         VectorField<TSScalar,SCOMPONENTS>    mImposed;

         /**
          * @brief Spectral vector total field
          */
         VectorField<TSScalar,SCOMPONENTS>    mTotal;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::total() const
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rTotal()
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::imposed() const
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS>& ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rImposed()
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

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void ImposedVectorVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::initialiseZeros()
   {
      // initialise the perturbation field to zero
      VectorVariable::initialiseZeros();

      // initialise vector imposed field to zero
      this->rImposed().initialiseZeros();

      // initialise vector total field to zero
      this->mTotal.initialiseZeros();
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
