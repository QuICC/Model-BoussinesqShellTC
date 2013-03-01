/** \file StorageDispatcher.hpp
 *  \brief Simple storage dispatcher for the communicator storage
 *
 *  \mhdBug Needs test
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

namespace GeoMHDiSCC {

namespace Parallel {

   template <Dimensions::Transform::Id TID> class StorageDispatcher;

   template <> class StorageDispatcher<Dimensions::Transform::TRA1D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProvider<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };
   
   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProvider<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA1D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage1D;
   }

   template <> class StorageDispatcher<Dimensions::Transform::TRA2D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProvider<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProvider<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA2D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage2D;
   }

   template <> class StorageDispatcher<Dimensions::Transform::TRA3D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static StoragePairProvider<typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> StoragePairProvider<typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  StorageDispatcher<Dimensions::Transform::TRA3D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.mStorage3D;
   }
}
}

#endif // STORAGEDISPATCHER_HPP
