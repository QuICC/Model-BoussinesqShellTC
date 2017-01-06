/** 
 * @file VectorVariable.hpp
 * @brief Implementation of vector variable
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

namespace QuICC {

namespace Datatypes {

   /**
    * @brief Implementation of vector variable
    */
   template <typename TSScalar, typename TPScalar> class VectorVariable: public VectorPhysicalVariable<TPScalar>
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
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  perturbation() const;

         /**
          * @brief Set spectral vector field (perturbation part)
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rPerturbation();

         /**
          * @brief Get spectral vector field (total field)
          */
         const VectorField<TSScalar,FieldComponents::Spectral::Id>&  total() const;

         /**
          * @brief Set spectral vector field (total field)
          */
         VectorField<TSScalar,FieldComponents::Spectral::Id>&  rTotal();

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
          * @brief Spectral vector of the field
          */
         SharedPtrMacro<VectorField<TSScalar,FieldComponents::Spectral::Id> > mspPerturbation;

      private:
   };

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::perturbation() const
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::rPerturbation()
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline const VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::total() const
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> inline VectorField<TSScalar,FieldComponents::Spectral::Id>& VectorVariable<TSScalar,TPScalar>::rTotal()
   {
      return *this->mspPerturbation;
   }

   template <typename TSScalar, typename TPScalar> VectorVariable<TSScalar,TPScalar>::VectorVariable(SharedResolution spRes)
      : VectorPhysicalVariable<TPScalar>(spRes)
   {
   }

   template <typename TSScalar, typename TPScalar> VectorVariable<TSScalar,TPScalar>::~VectorVariable()
   {
   }

   template <typename TSScalar, typename TPScalar> void VectorVariable<TSScalar,TPScalar>::setZeros()
   {
      // initialise the physical components to zero
      VectorPhysicalVariable<TPScalar>::setZeros();

      // initialise vector field to zero
      this->rPerturbation().setZeros();
   }

   template <typename TSScalar, typename TPScalar> void VectorVariable<TSScalar,TPScalar>::initSpectral(const std::vector<FieldComponents::Spectral::Id>& comps)
   {
      // Safety assert
      assert(! this->mspPerturbation);

      std::map<FieldComponents::Spectral::Id,bool> map;
      for(unsigned int i = 0; i < comps.size(); i++)
      {
         map.insert(std::make_pair(comps.at(i), true));
      }

      this->mspPerturbation = SharedPtrMacro<VectorField<TSScalar,FieldComponents::Spectral::Id> >(new VectorField<TSScalar,FieldComponents::Spectral::Id>(this->spRes()->spSpectralSetup(), map));
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TSScalar, typename TPScalar> MHDFloat VectorVariable<TSScalar,TPScalar>::requiredStorage() const
   {
      MHDFloat mem = 0.0;

      mem += VectorPhysicalVariable<TPScalar>::requiredStorage();

      mem += this->mspPerturbation->requiredStorage();

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // VECTORVARIABLE_HPP
