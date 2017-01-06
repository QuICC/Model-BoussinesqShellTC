/** 
 * @file SymmetricTensorField.hpp
 * @brief Implementation of a generic symmetric tensor field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SYMMETRICTENSORFIELD_HPP
#define SYMMETRICTENSORFIELD_HPP

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
    * @brief Implementation of a generic symmetric tensor field
    */
   template <typename TScalar, typename TType> class SymmetricTensorField
   {
      public:
         /**
          * @brief Constructor with simplified interface
          *
          * @param spSetup Shared setup object for the scalar fields
          */
         SymmetricTensorField(typename TScalar::SharedSetupType spSetup, const typename std::map<std::pair<TType,TType>,bool>& comps);

         /**
          * @brief Destructor
          */
         virtual ~SymmetricTensorField();

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

      #ifdef QUICC_STORAGEPROFILE
         /**
          * @brief Get the memory requirements
          */
         MHDFloat requiredStorage() const;
     #endif // QUICC_STORAGEPROFILE
         
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

   template <typename TScalar, typename TType> inline const TScalar& SymmetricTensorField<TScalar,TType>::comp(const TType id1, const TType id2) const
   {
      // Assert that index is valid
      assert(this->mComponents.count(std::make_pair(id1, id2)) == 1 || this->mComponents.count(std::make_pair(id2, id1)) == 1);

      if(this->mComponents.count(std::make_pair(id1, id2)) == 1)
      {
         return this->mComponents.find(std::make_pair(id1,id2))->second;
      } else
      {
         return this->mComponents.find(std::make_pair(id2,id1))->second;
      }
   }

   template <typename TScalar, typename TType> inline TScalar& SymmetricTensorField<TScalar,TType>::rComp(const TType id1, const TType id2)
   {
      // Assert that index is valid
      assert(this->mComponents.count(std::make_pair(id1,id2)) == 1 || this->mComponents.count(std::make_pair(id2,id1)) == 1);

      if(this->mComponents.count(std::make_pair(id1,id2)) == 1)
      {
         return this->mComponents.find(std::make_pair(id1,id2))->second;
      } else
      {
         return this->mComponents.find(std::make_pair(id2,id1))->second;
      }
   }
   
   template <typename TScalar, typename TType> inline const std::map<std::pair<TType,TType>,TScalar>& SymmetricTensorField<TScalar,TType>::data() const
   {
      return this->mComponents;
   }
   
   template <typename TScalar, typename TType> inline std::map<std::pair<TType,TType>,TScalar>& SymmetricTensorField<TScalar,TType>::rData()
   {
      return this->mComponents;
   }

   template <typename TScalar, typename TType> SymmetricTensorField<TScalar,TType>::SymmetricTensorField(typename TScalar::SharedSetupType spSetup, const typename std::map<std::pair<TType,TType>,bool>& comps)
   {
      // Initialise the components
      typename std::map<std::pair<TType,TType>,bool>::const_iterator   it;
      for(it = comps.begin(); it != comps.end(); ++it)
      {
         // Symmetric tensor so only sub or super diagonals
         if(it->second && this->mComponents.count(std::make_pair(it->first.second, it->first.first)) == 0)
         {
            this->mComponents.insert(std::make_pair(it->first, TScalar(spSetup)));
         }
      }
   }

   template <typename TScalar, typename TType> SymmetricTensorField<TScalar, TType>::~SymmetricTensorField()
   {
   }

   template <typename TScalar, typename TType> void SymmetricTensorField<TScalar,TType>::setZeros()
   {
      // Initialise the components
      typename std::map<std::pair<TType,TType>,TScalar>::iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         it->second.setZeros();
      }
   }

#ifdef QUICC_STORAGEPROFILE
   template <typename TScalar, typename TType> MHDFloat SymmetricTensorField<TScalar,TType>::requiredStorage() const
   {
      MHDFloat mem = 0.0;
      typename std::map<std::pair<TType,TType>,TScalar>::const_iterator   it;
      for(it = this->mComponents.begin(); it != this->mComponents.end(); it++)
      {
         mem += it->second.requiredStorage();
      }

      return mem;
   }
#endif // QUICC_STORAGEPROFILE

}
}

#endif // SYMMETRICTENSORFIELD_HPP
