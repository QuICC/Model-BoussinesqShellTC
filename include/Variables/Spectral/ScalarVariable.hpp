/** 
 * @file ScalarVariable.hpp
 * @brief Implementation of scalar field variable
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SCALARVARIABLE_HPP
#define SCALARVARIABLE_HPP

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
#include "Variables/Physical/ScalarPhysicalVariable.hpp"

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of scalar field variable
    */
   template <typename TSScalar, typename TPScalar> class ScalarVariable: public ScalarPhysicalVariable<TPScalar>
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         ScalarVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~ScalarVariable();

         /**
          * @brief Get scalar field (perturbation part)
          */
         const TSScalar&  perturbation() const;

         /**
          * @brief Set scalar field (perturbation part)
          */
         TSScalar&  rPerturbation();

         /**
          * @brief Get scalar field (total field)
          */
         const TSScalar&  total() const;

         /**
          * @brief Set scalar field (total field)
          */
         TSScalar&  rTotal();

         /**
          * @brief initialise to zeros
          */
         void setZeros();

         /**
          * @brief Initialise the spectral values storage
          */
         void initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps);

     #ifdef QUICC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Spectral scalar field
          */
         TSScalar    mPerturbation;

      private:
   };

   template <typename TSScalar, typename TPScalar> inline const TSScalar& ScalarVariable<TSScalar,TPScalar>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline TSScalar& ScalarVariable<TSScalar,TPScalar>::rPerturbation()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline const TSScalar& ScalarVariable<TSScalar,TPScalar>::total() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline TSScalar& ScalarVariable<TSScalar,TPScalar>::rTotal()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, typename TPScalar> ScalarVariable<TSScalar,TPScalar>::ScalarVariable(SharedResolution spRes)
      : ScalarPhysicalVariable<TPScalar>(spRes), mPerturbation(spRes->spSpectralSetup())
   {
   }

   template <typename TSScalar, typename TPScalar> ScalarVariable<TSScalar,TPScalar>::~ScalarVariable()
   {
   }

   template <typename TSScalar, typename TPScalar> void ScalarVariable<TSScalar,TPScalar>::setZeros()
   {
      // initialise the physical components to zero
      ScalarPhysicalVariable<TPScalar>::setZeros();

      // initialise scalar field to zero
      this->rPerturbation().setZeros();
   }

   template <typename TSScalar, typename TPScalar> void ScalarVariable<TSScalar,TPScalar>::initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps)
   {
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TSScalar, typename TPScalar> MHDFloat ScalarVariable<TSScalar,TPScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += ScalarPhysicalVariable<TPScalar>::requiredStorage();

      mem += this->mPerturbation.requiredStorage();

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // SCALARVARIABLE_HPP
