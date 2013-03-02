/** \file ScalarVariable.hpp
 *  \brief Implementation of scalar field variable
 *
 *  \mhdBug Needs test
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

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * \brief Implementation of scalar field variable
    */
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> class ScalarVariable: public ScalarPhysicalVariable<TPScalar, PCOMPONENTS>
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
         void initialiseZeros();

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:
         /**
          * @brief Spectral scalar field
          */
         TSScalar    mPerturbation;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::perturbation() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::rPerturbation()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::total() const
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::rTotal()
   {
      return this->mPerturbation;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::ScalarVariable(SharedResolution spRes)
      : ScalarPhysicalVariable<TPScalar,PCOMPONENTS>(spRes), mPerturbation(*spRes->spSpectralSetup())
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::~ScalarVariable()
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::initialiseZeros()
   {
      // initialise the physical components to zero
      ScalarPhysicalVariable<TPScalar,PCOMPONENTS>::initialiseZeros();

      // initialise scalar field to zero
      this->rPerturbation().initialiseZeros();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> MHDFloat ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += ScalarPhysicalVariable<TPScalar,PCOMPONENTS>::requiredStorage();

      mem += this->mPerturbation.requiredStorage();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // SCALARVARIABLE_HPP
