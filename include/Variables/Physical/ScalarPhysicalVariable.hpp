/** 
 * @file ScalarPhysicalVariable.hpp
 * @brief Base of the implementation of the physical components of a scalar variable
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
   template <typename TScalar> class ScalarPhysicalVariable : public VariableBase
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
         const VectorField<TScalar,FieldComponents::Physical::Id>&   grad() const;

         /**
          * @brief Set the physical gradient values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rGrad();

         /**
          * @brief Initialise to zero
          */
         void setZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical gradient storage
          */
         void initPhysicalGradient(const FieldComponents::Physical::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Check if variable has physical data setup
          */
         bool hasPhys() const;

         /**
          * @brief Check if variable has gradient data setup
          */
         bool hasGrad() const;

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
         SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > mspGrad;
   };

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasPhys() const
   {
      return this->mspPhys;
   }

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasGrad() const
   {
      return this->mspGrad;
   }

   template <typename TScalar> inline const TScalar&  ScalarPhysicalVariable<TScalar>::phys() const
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline TScalar&  ScalarPhysicalVariable<TScalar>::rPhys()
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::grad() const
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::rGrad()
   {
      // Safety assertion
      assert(this->mspGrad);

      return *this->mspGrad;
   }

   template <typename TScalar> ScalarPhysicalVariable<TScalar>::ScalarPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar> ScalarPhysicalVariable<TScalar>::~ScalarPhysicalVariable()
   {
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::setZeros()
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

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspPhys);

      this->mspPhys = SharedPtrMacro<TScalar>(new TScalar(this->spRes()->spPhysicalSetup()));
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysicalGradient(const FieldComponents::Physical::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspGrad);
      assert(id == FieldComponents::Physical::NOTUSED);

      this->mspGrad = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar> MHDFloat ScalarPhysicalVariable<TScalar>::requiredStorage() const
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
