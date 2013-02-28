/** \file VectorVariable.hpp
 *  \brief Implementation of vector variable
 *
 *  \mhdBug Needs test
 */

#ifndef VECTORVARIABLE_HPP
#define VECTORVARIABLE_HPP

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
#include "VectorFields/VectorField.hpp"
#include "Variables/Physical/VectorPhysicalVariable.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * \brief Implementation of vector variable
    */
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> class VectorVariable: public VectorPhysicalVariable<TPScalar, PCOMPONENTS>
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         VectorVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~VectorVariable();

         /**
          * @brief Get spectral vector field (perturbation part)
          */
         const VectorField<TSScalar,SCOMPONENTS>&  perturbation() const;

         /**
          * @brief Set spectral vector field (perturbation part)
          */
         VectorField<TSScalar,SCOMPONENTS>&  rPerturbation();

         /**
          * @brief Get spectral vector field (total field)
          */
         const VectorField<TSScalar,SCOMPONENTS>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,SCOMPONENTS>&  rTotal();

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
          * @brief Spectral vector of the field
          */
         VectorField<TSScalar,SCOMPONENTS>   mPerturbation;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS>& VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS>& VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::rPerturbation()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const VectorField<TSScalar,SCOMPONENTS>& VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::total() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline VectorField<TSScalar,SCOMPONENTS>& VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::rTotal()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::VectorVariable(SharedResolution spRes)
      : VectorPhysicalVariable<TPScalar,PCOMPONENTS>(spRes), mPerturbation(*spRes->spBwdSetup())
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::~VectorVariable()
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::initialiseZeros()
   {
      // initialise the physical components to zero
      VectorPhysicalVariable<TPScalar,PCOMPONENTS>::initialiseZeros();

      // initialise vector field to zero
      this->rPerturbation().initialiseZeros();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> MHDFloat VectorVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += VectorPhysicalVariable<TPScalar,PCOMPONENTS>::requiredStorage();

      mem += this->mPerturbation.requiredStorage();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // VECTORVARIABLE_HPP
