/** \file CommunicatorStorage.hpp
 *  \brief Simple implementation of the storage providers for the communicators
 *
 *  \mhdBug Needs test
 */

#ifndef COMMUNICATORSTORAGE_HPP
#define COMMUNICATORSTORAGE_HPP

// Configuration includes
//
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <vector>

// External includes
//

// Project includes
//
#include "StorageProviders/StoragePairProvider.hpp"
#include "Communicators/CommunicatorBase.hpp"
#include "TypeSelectors/IndexConverterSelector.hpp"
#include "Communicators/Converters/IConverter.hpp"
#include "Communicators/Converters/SerialConverter.hpp"
#ifdef GEOMHDISCC_MPI
   #include "Communicators/Converters/MpiConverter.hpp"
#endif // GEOMHDISCC_MPI

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Simple implementation of the storage providers for the communicators
    */ 
   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> class CommunicatorStorage;

   /**
    * @brief Specialisation of the implementation of the storage providers for the communicators
    */ 
   template <template <Dimensions::Transform::Id> class TTypes> class CommunicatorStorage<Dimensions::ONED, TTypes>: public CommunicatorBase
   {
      public:
         /// Typedef for the forward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::FwdType Fwd1DType;

         /// Typedef for the backward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::BwdType Bwd1DType;

         /**
         * @brief Constructor
         */
         CommunicatorStorage() {};

         /**
         * @brief Destructor
         */
         virtual ~CommunicatorStorage() {};

         /**
          * @brief Get/Set storage provider
          */
         template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& storage();

         /**
          * @brief Get/Set data converter
          */
         template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType, typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& converter();
         
      protected:

      private:
         /**
          * @brief Storage provider for 1D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>  mStorage1D; 
   };

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::ONED,TTypes>::storage()
   {
      Debug::StaticAssert< static_cast<int>(TID) <= static_cast<int>(Dimensions::ONED) >();

      return this->mStorage1D;
   }

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType,typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::ONED,TTypes>::converter()
   {
      Debug::StaticAssert< false >();

      return -1;
   }

   /**
    * @brief Specialisation of the implementation of the storage providers for 2D the communicators
    */ 
   template <template <Dimensions::Transform::Id> class TTypes> class CommunicatorStorage<Dimensions::TWOD, TTypes>: public CommunicatorBase
   {
      public:
         /// Typedef for the forward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::FwdType Fwd1DType;

         /// Typedef for the backward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::BwdType Bwd1DType;

         /// Typedef for the forward data type in 2D
         typedef typename TTypes<Dimensions::Transform::TRA2D>::FwdType Fwd2DType;

         /// Typedef for the backward data type in 2D
         typedef typename TTypes<Dimensions::Transform::TRA2D>::BwdType Bwd2DType;

         /**
         * @brief Constructor
         */
         CommunicatorStorage() {};

         /**
         * @brief Destructor
         */
         virtual ~CommunicatorStorage() {};

         /**
          * @brief Get/Set storage provider
          */
         template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& storage();

         /**
          * @brief Get/Set data converter
          */
         template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType, typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& converter();
         
      protected:
         /**
          * @brief Converter between first and second transform
          */
         SharedPtrMacro<IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType> > mspConverter2D;

      private:
         /**
          * @brief Storage provider for 1D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>  mStorage1D; 

         /**
          * @brief Storage provider for 2D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>  mStorage2D; 
   };

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::TWOD,TTypes>::storage()
   {
      Debug::StaticAssert< static_cast<int>(TID) <= static_cast<int>(Dimensions::TWOD) >();

      switch(TID)
      {
         case(Dimensions::Transform::TRA1D):
            return this->mStorage1D;
         case(Dimensions::Transform::TRA2D):
            return this->mStorage2D;
      }
   }

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType, typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::TWOD,TTypes>::converter()
   {
      Debug::StaticAssert< TID == Dimensions::Transform::TRA2D >();

      return *this->mspConverters2D;
   }

   /**
    * @brief Specialisation of the implementation of the storage providers for the 3D communicators
    */ 
   template <template <Dimensions::Transform::Id> class TTypes> class CommunicatorStorage<Dimensions::THREED, TTypes>: public CommunicatorBase
   {
      public:
         /// Typedef for the forward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::FwdType Fwd1DType;

         /// Typedef for the backward data type in 1D
         typedef typename TTypes<Dimensions::Transform::TRA1D>::BwdType Bwd1DType;

         /// Typedef for the forward data type in 2D
         typedef typename TTypes<Dimensions::Transform::TRA2D>::FwdType Fwd2DType;

         /// Typedef for the backward data type in 2D
         typedef typename TTypes<Dimensions::Transform::TRA2D>::BwdType Bwd2DType;

         /// Typedef for the forward data type in 3D
         typedef typename TTypes<Dimensions::Transform::TRA3D>::FwdType Fwd3DType;

         /// Typedef for the backward data type in 3D
         typedef typename TTypes<Dimensions::Transform::TRA3D>::BwdType Bwd3DType;

         /**
         * @brief Constructor
         */
         CommunicatorStorage() {};

         /**
         * @brief Destructor
         */
         virtual ~CommunicatorStorage() {};

         /**
          * @brief Get/Set storage provider
          */
         template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& storage();

         /**
          * @brief Get/Set data converter
          */
         template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType, typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& converter();
         
      protected:
         /**
          * @brief Converter between first and second transform
          */
         SharedPtrMacro<IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType> > mspConverter2D;

         /**
          * @brief Converter between second and third transform
          */
         SharedPtrMacro<IConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType> > mspConverter3D;

      private:
         /**
          * @brief Storage provider for 1D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>  mStorage1D; 

         /**
          * @brief Storage provider for 2D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>  mStorage2D; 

         /**
          * @brief Storage provider for 3D
          */
         StoragePairProvider<typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>  mStorage3D; 
   };

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> StoragePairProvider<typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::THREED,TTypes>::storage()
   {
      Debug::StaticAssert< static_cast<int>(TID) <= static_cast<int>(Dimensions::THREED) >();

      switch(TID)
      {
         case(Dimensions::Transform::TRA1D):
            return this->mStorage1D;
         case(Dimensions::Transform::TRA2D):
            return this->mStorage2D;
         case(Dimensions::Transform::TRA3D):
            return this->mStorage3D;
      }
   }

   template <template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> IConverter<typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::FwdType,typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(TID)-1)>::BwdType, typename TTypes<TID>::FwdType,typename TTypes<TID>::BwdType>& CommunicatorStorage<Dimensions::THREED,TTypes>::converter()
   {
      Debug::StaticAssert< TID != Dimensions::Transform::TRA1D >();

      switch(TID)
      {
         case(Dimensions::Transform::TRA2D):
            return *this->mspConverters2D;
         case(Dimensions::Transform::TRA3D):
            return *this->mspConverters3D;
      }
   }

}
}

#endif // COMMUNICATORSTORAGE_HPP
