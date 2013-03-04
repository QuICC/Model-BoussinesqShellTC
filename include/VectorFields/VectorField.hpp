/** \file VectorField.hpp
 *  \brief Implementation of a generic vector field
 *
 *  \mhdBug Needs test
 */

#ifndef VECTORFIELD_HPP
#define VECTORFIELD_HPP

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldComponents.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Implementation of a generic vector field
    */
   template <typename TScalar, int COMPONENTS, typename TType> class VectorField
   {
      public:
         /**
          * @brief Constructor with simplified interface
          *
          * @param spSetup Shared setup object for the scalar fields
          */
         VectorField(typename TScalar::SharedSetupType spSetup);

         /**
          * @brief Destructor
          */
         virtual ~VectorField();

         /**
          * @brief Get field component
          *
          * @param i Index of the component
          */
         const TScalar& comp(const TType id) const;

         /**
          * @brief Set field component
          *
          * @param i Index of the component
          */
         TScalar& rComp(const TType id);

         /**
          * @brief Get field components
          */
         const std::vector<TScalar>& data() const;

         /**
          * @brief Set field components to zero
          */
         void setZeros();

      #ifdef GEOMHDISCC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // GEOMHDISCC_STORAGEPROFILE

         /**
          * @brief Set internal storage field data
          *
          * \warning This routine should only be used in exceptional cases. Use setData, addData, subData when you can!
          */
         std::vector<TScalar>& rData();
         
      protected:
         /**
          * @brief Storage for the field components
          */
         std::vector<TScalar> mComponents;

      private:
   };

   template <typename TScalar, int COMPONENTS, typename TType> inline const TScalar& VectorField<TScalar,COMPONENTS,TType>::comp(const TType id) const
   {
      // Assert that index is valid
      assert(this->mComponents.size() > static_cast<size_t>(id));

      return this->mComponents.at(static_cast<int>(id));
   }

   template <typename TScalar, int COMPONENTS, typename TType> inline TScalar& VectorField<TScalar,COMPONENTS,TType>::rComp(const TType id)
   {
      // Assert that index is valid
      assert(this->mComponents.size() > static_cast<size_t>(id));

      return this->mComponents.at(static_cast<int>(id));
   }

   template <typename TScalar, int COMPONENTS, typename TType> VectorField<TScalar, COMPONENTS, TType>::VectorField(typename TScalar::SharedSetupType spSetup)
   {
      // Initialise the components
      for(int i = 0; i < COMPONENTS; i++)
      {
         this->mComponents.push_back(TScalar(spSetup));
      }
   }

   template <typename TScalar, int COMPONENTS, typename TType> VectorField<TScalar, COMPONENTS, TType>::~VectorField()
   {
   }

   template <typename TScalar, int COMPONENTS, typename TType> inline const std::vector<TScalar>& VectorField<TScalar,COMPONENTS,TType>::data() const
   {
      return this->mComponents;
   }

   template <typename TScalar, int COMPONENTS, typename TType> inline std::vector<TScalar>& VectorField<TScalar,COMPONENTS,TType>::rData()
   {
      return this->mComponents;
   }

   template <typename TScalar, int COMPONENTS, typename TType> void VectorField<TScalar, COMPONENTS, TType>::setZeros()
   {
      // Initialise the components
      typename std::vector<TScalar>::iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         it->setZeros();
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar, int COMPONENTS, typename TType> MHDFloat VectorField<TScalar, COMPONENTS, TType>::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      typename std::vector<TScalar>::const_iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         mem += it->requiredStorage();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // VECTORFIELD_HPP
