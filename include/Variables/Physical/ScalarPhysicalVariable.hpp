/** \file ScalarPhysicalVariable.hpp
 *  \brief Base of the implementation of the physical components of a scalar variable
 */

#ifndef SCALARPHYSICALVARIABLE_HPP
#define SCALARPHYSICALVARIABLE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"
#include "StorageProfiler/StorageProfilerMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"
#include "VectorFields/VectorField.hpp"
#include "Variables/VariableBase.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Base of the implementation of the physical components of a scalar variable
    */
   template <typename TScalar, int COMPONENTS> class ScalarPhysicalVariable : public VariableBase
   {
      public:
         /**
         * @brief Constructor
         *
         * There is only very little done directly. The different fields have to be
         * initialised after construction.
         *
         * @param spRes Resolution information
         */
         ScalarPhysicalVariable(SharedResolution spRes);

         /**
         * @brief Simple empty destructor
         */
         virtual ~ScalarPhysicalVariable();

         /**
          * @brief Get the physical field values
          */
         const TScalar&   phys() const;

         /**
          * @brief Set the physical field values
          */
         TScalar&   rPhys();

         /**
          * @brief Get the physical gradient values
          */
         const VectorField<TScalar,COMPONENTS>&   grad() const;

         /**
          * @brief Set the physical gradient values
          */
         VectorField<TScalar,COMPONENTS>&   rGrad();

         /**
          * @brief Initialise to zero
          */
         void initialiseZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initialisePhysical();

         /**
          * @brief Initialise the physical gradient storage
          */
         void initialisePhysicalDiff();

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Smart pointer for the physical field values
          */
         SharedPtrMacro<TScalar> mspPhys;

         /**
          * @brief Smart pointer for the physical gradient values
          */
         SharedPtrMacro<VectorField<TScalar,COMPONENTS> > mspGrad;
   };

   template <typename TScalar, int COMPONENTS> inline const TScalar&  ScalarPhysicalVariable<TScalar,COMPONENTS>::phys() const
   {
      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline TScalar&  ScalarPhysicalVariable<TScalar,COMPONENTS>::rPhys()
   {
      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline const VectorField<TScalar,COMPONENTS>&  ScalarPhysicalVariable<TScalar,COMPONENTS>::grad() const
   {
      return *this->mspGrad;
   }

   template <typename TScalar, int COMPONENTS> inline VectorField<TScalar,COMPONENTS>&  ScalarPhysicalVariable<TScalar,COMPONENTS>::rGrad()
   {
      return *this->mspGrad;
   }

   template <typename TScalar, int COMPONENTS> ScalarPhysicalVariable<TScalar,COMPONENTS>::ScalarPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar, int COMPONENTS> ScalarPhysicalVariable<TScalar,COMPONENTS>::~ScalarPhysicalVariable()
   {
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::initialiseZeros()
   {
      // Initialise physical values to zero if required
      if(this->mspPhys)
      {
         this->rPhys().initialiseZeros();
      }

      // Initialise gradient values to zero if required
      if(this->mspGrad)
      {
         this->rGrad().initialiseZeros();
      }
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::initialisePhysical()
   {
      this->mspPhys = SharedPtrMacro<TScalar>(new TScalar(*this->spRes()->spFwdSetup()));
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::initialisePhysicalDiff()
   {
      this->mspGrad = SharedPtrMacro<VectorField<TScalar,COMPONENTS> >(new VectorField<TScalar,COMPONENTS>(*this->spRes()->spFwdSetup()));
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar, int COMPONENTS> MHDFloat ScalarPhysicalVariable<TScalar,COMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Physical storage
      if(this->mspPhys)
      {
         mem += this->phys().requiredStorage();
      }

      // Physical gradient storage
      if(this->mspGrad)
      {
         mem += this->grad().requiredStorage();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE
}

#endif // SCALARPHYSICALVARIABLE_HPP
