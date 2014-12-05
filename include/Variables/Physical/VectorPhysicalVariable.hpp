/** 
 * @file VectorPhysicalVariable.hpp
 * @brief Base of the implementation of the physical components of a vector variable
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VECTORPHYSICALVARIABLE_HPP
#define VECTORPHYSICALVARIABLE_HPP

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
    * @brief Base of the implementation of the physical components of a vector variable
    */
   template <typename TScalar> class VectorPhysicalVariable : public VariableBase
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
         VectorPhysicalVariable(SharedResolution spRes);

         /**
         * @brief Simple empty destructor
         */
         virtual ~VectorPhysicalVariable();

         /**
          * @brief Get the physical field values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   phys() const;

         /**
          * @brief Set the physical field values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rPhys();

         /**
          * @brief Get the physical vector gradient values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   grad(FieldComponents::Physical::Id id) const;

         /**
          * @brief Set the physical vector gradient values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rGrad(FieldComponents::Physical::Id id);

         /**
          * @brief Get the physical curl values
          */
         const VectorField<TScalar,FieldComponents::Physical::Id>&   curl() const;

         /**
          * @brief Set the physical curl values
          */
         VectorField<TScalar,FieldComponents::Physical::Id>&   rCurl();

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
          * @brief Initialise the physical curl storage
          */
         void initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps);

         /**
          * @brief Check if variable has physical data setup
          */
         bool hasPhys() const;

         /**
          * @brief Check if variable has vector gradient data setup
          */
         bool hasGrad() const;

         /**
          * @brief Check if variable has curl data setup
          */
         bool hasCurl() const;

     #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
         * @brief Get the memory requirements
         */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE
         
      protected:

      private:
         /**
          * @brief Smart pointer for the physical space field values
          */
         SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > mspPhys;

         /**
          * @brief Smart pointer for the physical vector gradient values
          */
         std::map<FieldComponents::Physical::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > > mVGrad;

         /**
          * @brief Smart pointer for the physical curl values
          */
         SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > mspCurl;
   };

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasPhys() const
   {
      return this->mspPhys;
   }

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasGrad() const
   {
      return this->mVGrad.size();
   }

   template <typename TScalar> inline bool VectorPhysicalVariable<TScalar>::hasCurl() const
   {
      return this->mspCurl;
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::phys() const
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rPhys()
   {
      // Safety assertion
      assert(this->mspPhys);

      return *this->mspPhys;
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::grad(FieldComponents::Physical::Id id) const
   {
      // Safety assertion
      assert(this->mVGrad.count(id));

      return *(this->mVGrad.find(id)->second);
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rGrad(FieldComponents::Physical::Id id)
   {
      // Safety assertion
      assert(this->mVGrad.count(id));

      return *(this->mVGrad.find(id)->second);
   }

   template <typename TScalar> inline const VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::curl() const
   {
      // Safety assertion
      assert(this->mspCurl);

      return *this->mspCurl;
   }

   template <typename TScalar> inline VectorField<TScalar,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar>::rCurl()
   {
      // Safety assertion
      assert(this->mspCurl);

      return *this->mspCurl;
   }

   template <typename TScalar> VectorPhysicalVariable<TScalar>::VectorPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar> VectorPhysicalVariable<TScalar>::~VectorPhysicalVariable()
   {
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::setZeros()
   {
      // Initialise physical values to zero if required
      if(this->mspPhys)
      {
         this->rPhys().setZeros();
      }

      // Initialise vector gradient values to zero if required
      if(this->mVGrad.size() > 0)
      {
         typename std::map<FieldComponents::Physical::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > >::iterator it;
         for(it = this->mVGrad.begin(); it != this->mVGrad.end(); ++it)
         {
            it->second->setZeros();
         }
      }

      // Initialise curl values to zero if required
      if(this->mspCurl)
      {
         this->rCurl().setZeros();
      }
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysical(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspPhys);

      this->mspPhys = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysicalGradient(const FieldComponents::Physical::Id id, const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(this->mVGrad.count(id) == 0);

      // Create shared pointer
      SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > spGrad = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));

      // Insert into map
      this->mVGrad.insert(std::make_pair(id, spGrad));
   }

   template <typename TScalar> void VectorPhysicalVariable<TScalar>::initPhysicalCurl(const std::map<FieldComponents::Physical::Id,bool>& comps)
   {
      // Safety assert
      assert(! this->mspCurl);

      this->mspCurl = SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> >(new VectorField<TScalar,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup(), comps));
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar> MHDFloat VectorPhysicalVariable<TScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Physical storage
      if(this->mspPhys)
      {
         mem += this->phys().requiredStorage();
      }

      // Physical vector gradient storage
      if(this->mVGrad.size() > 0)
      {
         std::map<FieldComponents::Physical::Id,SharedPtrMacro<VectorField<TScalar,FieldComponents::Physical::Id> > >::const_iterator it;
         for(it = this->mVGrad.begin(); it != this->mVGrad.end(); ++it)
         {
            mem += it->second->requiredStorage();
         }
      }

      // Physical curl storage
      if(this->mspCurl)
      {
         mem += this->curl().requiredStorage();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE
}
}

#endif // VECTORPHYSICALVARIABLE_HPP
