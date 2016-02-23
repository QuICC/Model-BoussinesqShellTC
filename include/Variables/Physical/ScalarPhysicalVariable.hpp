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
          * @brief Get the physical 2nd order gradient values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   grad2(FieldComponents::Spectral::Id id) const;

         /**
          * @brief Set the physical 2nd order gradient values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rGrad2(FieldComponents::Spectral::Id id);

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
         void initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Initialise the physical 2nd order gradient storage
          */
         void initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Check if variable has physical data setup
          */
         bool hasPhys() const;

         /**
          * @brief Check if variable has gradient data setup
          */
         bool hasGrad() const;

         /**
          * @brief Check if variable has 2nd order gradient data setup
          */
         bool hasGrad2() const;

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

         /**
          * @brief Smart pointer for the physical 2nd order gradient values
          */
         std::map<FieldComponents::Spectral::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > > mGrad2;
   };

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasPhys() const
   {
      return this->mspPhys;
   }

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasGrad() const
   {
      return this->mspGrad;
   }

   template <typename TScalar> inline bool ScalarPhysicalVariable<TScalar>::hasGrad2() const
   {
      return this->mGrad2.size();
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

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::grad2(FieldComponents::Spectral::Id id) const
   {
      // Safety assertion
      assert(this->mGrad2.count(id));

      return *(this->mGrad2.find(id)->second);
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  ScalarPhysicalVariable<TScalar>::rGrad2(FieldComponents::Spectral::Id id)
   {
      // Safety assertion
      assert(this->mGrad2.count(id));

      return *(this->mGrad2.find(id)->second);
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

      // Initialise vector gradient values to zero if required
      if(this->mGrad2.size() > 0)
      {
         typename std::map<FieldComponents::Spectral::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > >::iterator it;
         for(it = this->mGrad2.begin(); it != this->mGrad2.end(); ++it)
         {
            it->second->setZeros();
         }
      }
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspPhys);

      this->mspPhys = SharedPtrMacro<TScalar>(new TScalar(this->spRes()->spPhysicalSetup()));
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysicalGradient(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspGrad);
      assert(id == FieldComponents::Spectral::SCALAR);

      this->mspGrad = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));
   }

   template <typename TScalar> void ScalarPhysicalVariable<TScalar>::initPhysicalGradient2(const FieldComponents::Spectral::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(this->mGrad2.count(id) == 0);

      // Create shared pointer
      SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > spGrad2 = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));

      // Insert into map
      this->mGrad2.insert(std::make_pair(id, spGrad2));
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

      // Physical 2nd order gradient storage
      if(this->mGrad2.size() > 0)
      {
         std::map<FieldComponents::Spectral::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > >::const_iterator it;
         for(it = this->mGrad2.begin(); it != this->mGrad2.end(); ++it)
         {
            mem += it->second->requiredStorage();
         }
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE
}
}

#endif // SCALARPHYSICALVARIABLE_HPP
