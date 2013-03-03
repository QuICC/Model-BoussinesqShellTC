/** \file VectorPhysicalVariable.hpp
 *  \brief Base of the implementation of the physical components of a vector variable
 *
 *  \mhdBug Needs test
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
   template <typename TScalar, int COMPONENTS> class VectorPhysicalVariable : public VariableBase
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
         const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   phys() const;

         /**
          * @brief Set the physical field values
          */
         VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   rPhys();

         /**
          * @brief Get the physical curl values
          */
         const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   curl() const;

         /**
          * @brief Set the physical curl values
          */
         VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&   rCurl();

         /**
          * @brief Initialise to zero
          */
         void setZeros();

         /**
          * @brief Initialise the physical values storage
          */
         void initializePhysical();

         /**
          * @brief Initialise the physical curl storage
          */
         void initializePhysicalDiff();

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
         SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> > mspPhys;

         /**
          * @brief Smart pointer for the physical curl values
          */
         SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> > mspCurl;
   };

   template <typename TScalar, int COMPONENTS> inline const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar,COMPONENTS>::phys() const
   {
      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar,COMPONENTS>::rPhys()
   {
      return *this->mspPhys;
   }

   template <typename TScalar, int COMPONENTS> inline const VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar,COMPONENTS>::curl() const
   {
      return *this->mspCurl;
   }

   template <typename TScalar, int COMPONENTS> inline VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>&  VectorPhysicalVariable<TScalar,COMPONENTS>::rCurl()
   {
      return *this->mspCurl;
   }

   template <typename TScalar, int COMPONENTS> VectorPhysicalVariable<TScalar,COMPONENTS>::VectorPhysicalVariable(SharedResolution spRes)
      : VariableBase(spRes)
   {
   }

   template <typename TScalar, int COMPONENTS> VectorPhysicalVariable<TScalar,COMPONENTS>::~VectorPhysicalVariable()
   {
   }

   template <typename TScalar, int COMPONENTS> void VectorPhysicalVariable<TScalar,COMPONENTS>::setZeros()
   {
      // Initialise physical values to zero if required
      if(this->mspPhys)
      {
         this->rPhys().setZeros();
      }

      // Initialise curl values to zero if required
      if(this->mspCurl)
      {
         this->rCurl().setZeros();
      }
   }

   template <typename TScalar, int COMPONENTS> void VectorPhysicalVariable<TScalar,COMPONENTS>::initializePhysical()
   {
      this->mspPhys = SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> >(new VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup()));
   }

   template <typename TScalar, int COMPONENTS> void VectorPhysicalVariable<TScalar,COMPONENTS>::initializePhysicalDiff()
   {
      this->mspCurl = SharedPtrMacro<VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id> >(new VectorField<TScalar,COMPONENTS,FieldComponents::Physical::Id>(this->spRes()->spPhysicalSetup()));
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar, int COMPONENTS> MHDFloat VectorPhysicalVariable<TScalar,COMPONENTS>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      // Physical storage
      if(this->mspPhys)
      {
         mem += this->phys().requiredStorage();
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
