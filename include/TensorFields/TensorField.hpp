/** 
 * @file TensorField.hpp
 * @brief Implementation of a generic tensor field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TENSORFIELD_HPP
#define TENSORFIELD_HPP

// Configuration includes
//

// System includes
//
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"

namespace GeoMHDiSCC {

namespace Datatypes {

   /**
    * @brief Implementation of a generic tensor field
    */
   template <typename TScalar, typename TType> class TensorField
   {
      public:
         /**
          * @brief Constructor with simplified interface
          *
          * @param spSetup Shared setup object for the scalar fields
          */
         TensorField(typename TScalar::SharedSetupType spSetup, const typename std::map<std::pair<TType,TType>,bool>& comps);

         /**
          * @brief Destructor
          */
         virtual ~TensorField();

         /**
          * @brief Get field component
          *
          * @param id1 Index of the first dimension
          * @param id2 Index of the first dimension
          */
         const TScalar& comp(const TType id1, const TType id2) const;

         /**
          * @brief Set field component
          *
          * @param id1 Index of the first dimension
          * @param id2 Index of the first dimension
          */
         TScalar& rComp(const TType id1, const TType id2);

         /**
          * @brief Set field components to zero
          */
         void setZeros();

         /**
          * @brief Get field components
          */
         const std::map<std::pair<TType,TType>,TScalar>& data() const;

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
         std::map<std::pair<TType,TType>,TScalar>& rData();

      protected:
         /**
          * @brief Storage for the field components
          */
         std::map<std::pair<TType,TType>,TScalar> mComponents;

      private:
   };

   template <typename TScalar, typename TType> inline const TScalar& TensorField<TScalar,TType>::comp(const TType id1, const TType id2) const
   {
      // Assert that index is valid
      assert(this->mComponents.count(std::make_pair(id1, id2)) == 1);

      return this->mComponents.find(std::make_pair(id1,id2))->second;
   }

   template <typename TScalar, typename TType> inline TScalar& TensorField<TScalar,TType>::rComp(const TType id1, const TType id2)
   {
      // Assert that index is valid
      assert(this->mComponents.count(std::make_pair(id1,id2)) == 1);

      return this->mComponents.find(std::make_pair(id1,id2))->second;
   }
   
   template <typename TScalar, typename TType> inline const std::map<std::pair<TType,TType>,TScalar>& TensorField<TScalar,TType>::data() const
   {
      return this->mComponents;
   }
   
   template <typename TScalar, typename TType> inline std::map<std::pair<TType,TType>,TScalar>& TensorField<TScalar,TType>::rData()
   {
      return this->mComponents;
   }

   template <typename TScalar, typename TType> TensorField<TScalar,TType>::TensorField(typename TScalar::SharedSetupType spSetup, const typename std::map<std::pair<TType,TType>,bool>& comps)
   {
      // Initialise the components
      typename std::map<std::pair<TType,TType>,bool>::const_iterator   it;
      for(it = comps.begin(); it != comps.end(); ++it)
      {
         if(it->second)
         {
            this->mComponents.insert(std::make_pair(it->first, TScalar(spSetup)));
         }
      }
   }

   template <typename TScalar, typename TType> TensorField<TScalar, TType>::~TensorField()
   {
   }

   template <typename TScalar, typename TType> void TensorField<TScalar,TType>::setZeros()
   {
      // Initialise the components
      typename std::map<std::pair<TType,TType>,TScalar>::iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         it->second.setZeros();
      }
   }

#ifdef GEOMHDISCC_STORAGEPROFILE
   template <typename TScalar, typename TType> MHDFloat TensorField<TScalar,TType>::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      typename std::map<std::pair<TType,TType>,TScalar>::const_iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         mem += it->second.requiredStorage();
      }

      return mem;
   }
#endif // GEOMHDISCC_STORAGEPROFILE

}
}

#endif // TENSORFIELD_HPP
