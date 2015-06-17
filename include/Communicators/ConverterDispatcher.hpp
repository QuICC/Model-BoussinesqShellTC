/**
 * @file ConverterDispatcher.hpp
 * @brief Simple converter dispatcher for the communicators 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "TypeSelectors/IndexConverterSelector.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Templated class to provide access to a different converter for each dimension and different dimensionality
    */
   template <Dimensions::Transform::Id TID> class ConverterDispatcher;

   /**
    * @brief Specialised converter dispatcher for the first transform
    */
   template <> class ConverterDispatcher<Dimensions::Transform::TRA1D>
   {
      public:
         /**
          * @brief Get the correct converter
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>&  get(TComm<DIMENSION,TTypes>& comm);
   };
 
   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>&  ConverterDispatcher<Dimensions::Transform::TRA1D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      // This should never get called!
      Debug::StaticAssert<DIMENSION == -99>();

      return comm.mspConverter1D;
   }

   /**
    * @brief Specialised converter dispatcher for the second transform
    */
   template <> class ConverterDispatcher<Dimensions::Transform::TRA2D>
   {
      public:
         /**
          * @brief Get the correct converter
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType,typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>&  ConverterDispatcher<Dimensions::Transform::TRA2D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return *comm.mspConverter2D;
   }

   /**
    * @brief Specialised converter dispatcher for the third transform
    */
   template <> class ConverterDispatcher<Dimensions::Transform::TRA3D>
   {
      public:
         /**
          * @brief Get the correct converter
          */
         template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static IConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA3D>::Type>&  get(TComm<DIMENSION,TTypes>& comm);
   };

   template<Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> IConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType,typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType,typename TTypes<Dimensions::Transform::TRA3D>::BwdType, typename IndexConverterSelector<Dimensions::Transform::TRA3D>::Type>&  ConverterDispatcher<Dimensions::Transform::TRA3D>::get(TComm<DIMENSION,TTypes>& comm)
   {
      return *comm.mspConverter3D;
   }
}
}

#endif // CONVERTERDISPATCHER_HPP
