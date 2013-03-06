/** \file ScalarPhysicalVariable.hpp
 *  \brief Base of the implementation of the physical components of a scalar variable
 *
 *  \mhdBug Needs test
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

namespace Datatypes {

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
         const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   grad() const;

         /**
          * @brief Set the physical gradient values
          */
         VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   rGrad();

         /**
          * @brief Initialise to zero
          */
         void setZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initPhysical();

         /**
          * @brief Initialise the physical gradient storage
          */
         void initPhysicalDiff();

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
         SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> > mspGrad;
   };

   template <typename TScalar, int COMPONENTS> inline const TScalar&  ScalarPhysicalVariable<TScalar,COMPONENTS>::phys() const
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline TScalar&  ScalarPhysicalVariable<TScalar,COMPONENTS>::rPhys()
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar,COMPONENTS>::grad() const
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar, int COMPONENTS> inline VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar,COMPONENTS>::rGrad()
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar, int COMPONENTS> ScalarPhysicalVariable<TScalar,COMPONENTS>::ScalarPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar, int COMPONENTS> ScalarPhysicalVariable<TScalar,COMPONENTS>::~ScalarPhysicalVariable()
   {
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::setZeros()
   {
      // Initialise physical values to zero if required
      if(this->mspPhys)
      {
         this->rPhys().setZeros();
      }

      // Initialise gradient values to zero if required
      if(this->mspGrad)
      {
         this->rGrad().setZeros();
      }
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::initPhysical()
   {
      this->mspPhys = SharedPtrMacro<TScalar>(new TScalar(this->spRes()->spPhysicalSetup()));
   }

   template <typename TScalar, int COMPONENTS> void ScalarPhysicalVariable<TScalar,COMPONENTS>::initPhysicalDiff()
   {
      this->mspGrad = SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> >(new VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup()));
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
}

#endif // SCALARPHYSICALVARIABLE_HPP
