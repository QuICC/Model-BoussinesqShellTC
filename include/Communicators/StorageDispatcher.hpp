/**
 * @file StorageDispatcher.hpp
 * @brief Simple storage dispatcher for the communicator storage 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STORAGEDISPATCHER_HPP
#define STORAGEDISPATCHER_HPP

// Configuration includes
//
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Template class to provide different storage dispatcher for each transform
    */
   template <Dimensions::Transform::Id TID> class StorageDispatcher;

   /**
    * @brief Specialised storage dispatcher for the first transform
    */
   template <> class StorageDispatcher<Dimensions::Transform::TRA1D>
   {
      public:
         /**
          * @brief Get the storage provider corresponding the transform
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };
   
   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA1D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage1D;
   }

   /**
    * @brief Specialised storage dispatcher for the second transform
    */
   template <> class StorageDispatcher<Dimensions::Transform::TRA2D>
   {
      public:
         /**
          * @brief Get the storage provider corresponding the transform
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA2D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage2D;
   }

   /**
    * @brief Specialised storage dispatcher for the third transform
    */
   template <> class StorageDispatcher<Dimensions::Transform::TRA3D>
   {
      public:
         /**
          * @brief Get the storage provider corresponding the transform
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProviderMacro<typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA3D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage3D;
   }
}
}

#endif // STORAGEDISPATCHER_HPP
