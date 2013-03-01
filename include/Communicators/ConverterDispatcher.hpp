/** \file ConverterDispatcher.hpp
 *  \brief Simple converter dispatcher for the communicators
 *
 *  \mhdBug Needs test
 */

#ifndef CONVERTERDISPATCHER_HPP
#define CONVERTERDISPATCHER_HPP

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
#include "Communicators/Converters/IConverter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   template <Dimensions::Transform::Id TID> class ConverterDispatcher;

   template <> class ConverterDispatcher<Dimensions::Transform::TRA1D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };
   
   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  ConverterDispatcher<Dimensions::Transform::TRA1D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      Debug::StaticAssert<false>();
      return comm.mspConverter1D;
   }

   template <> class ConverterDispatcher<Dimensions::Transform::TRA2D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType>&  ConverterDispatcher<Dimensions::Transform::TRA2D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return *comm.mspConverter2D;
   }

   template <> class ConverterDispatcher<Dimensions::Transform::TRA3D>
   {
      public:
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType>&  ConverterDispatcher<Dimensions::Transform::TRA3D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return *comm.mspConverter3D;
   }
}
}

#endif // CONVERTERDISPATCHER_HPP
