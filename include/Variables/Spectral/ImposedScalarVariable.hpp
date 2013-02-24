/** \file ImposedScalarVariable.hpp
 *  \brief Implementation of scalar field variable with imposed component
 */

#ifndef IMPOSEDSCALARVARIABLE_HPP
#define IMPOSEDSCALARVARIABLE_HPP

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
#include "Variables/Spectral/ScalarVariable.hpp"

namespace GeoMHDiSCC {

   /**
    * \brief Implementation of scalar field variable with imposed component
    */
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> class ImposedScalarVariable: public ScalarVariable<TSScalar,SCOMPONENTS,TPScalar,PCOMPONENTS> 
   {
      public:
         /**
          * @brief Constructs the underlying physical and spectral fields
          *
          * @param spRes Resolution information
          */
         ImposedScalarVariable(SharedResolution spRes);

         /**
         * @brief Destructor
         */
         virtual ~ImposedScalarVariable();

         /**
          * @brief Get scalar field (total field)
          */
         const TSScalar&  total() const;

         /**
          * @brief Set scalar field (total field)
          */
         TSScalar&  rTotal();

         /**
          * @brief Get the imposed scalar field
          */
         const TSScalar&  imposed() const;

         /**
          * @brief Set the imposed scalar field
          */
         TSScalar&  rImposed();

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
          * @brief Imposed spectral scalar field
          */
         TSScalar    mImposed;

         /**
          * @brief Total scalar field (imposed + perturbation)
          */
         TSScalar    mTotal;

      private:
   };

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::total() const
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rTotal()
   {
      return this->mTotal;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline const TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::imposed() const
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> inline TSScalar& ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::rImposed()
   {
      return this->mImposed;
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedScalarVariabl<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>e::ImposedScalarVariable(SharedResolution spRes)
      : ScalarVariable(spRes), mImposed(*spRes->spBwdSetup()), mTotal(*spRes->spBwdSetup())
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> ImposedScalarVariabl<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>e::~ImposedScalarVariable()
   {
   }

   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> void ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::initialiseZeros()
   {
      // initialise the pertubation component to zero
      ScalarVariable::initialiseZeros();

      // initialise imposed scalar field to zero
      this->rImposed().initialiseZeros();

      // initialise total scalar field to zero
      this->mTotal.initialiseZeros();
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TSScalar, int SCOMPONENTS, typename TPScalar, int PCOMPONENTS> MHDFloat ImposedScalarVariable<TSScalar,SCOMPONENTS, TPScalar, PCOMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += ScalarVariable::requiredStorage();

      mem += this->mImposed.requiredStorage();

      mem += this->mTotal.requiredStorage();

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}

#endif // IMPOSEDSCALARVARIABLE_HPP
