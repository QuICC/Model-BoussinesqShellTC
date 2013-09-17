/**
 * @file EndDispatcher.hpp
 * @brief Simple converter dispatcher for the communicators 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ENDDISPATCHER_HPP
#define ENDDISPATCHER_HPP

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

   template <bool, Dimensions::Transform::Id TID> class EndDispatcher;

   template <Dimensions::Transform::Id TID> class EndDispatcher<true,TID>
   {
      public:
         /**
          * @brief Receive forward the data
          */
         template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static typename TTypes<TID>::FwdType& receiveForward(TComm<DIMENSION,TTypes>& comm);

         /**
          * @brief transfer forward data
          */
         template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static void transferForward(TComm<DIMENSION,TTypes>& comm, typename TTypes<TID>::FwdType& rData);

   };

   template <Dimensions::Transform::Id TID> template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> typename TTypes<TID>::FwdType& EndDispatcher<true,TID>::receiveForward(TComm<DIMENSION,TTypes>& comm)
   {
      return comm.template storage<TID>().recoverFwd();
   }

   template <Dimensions::Transform::Id TID> template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> void EndDispatcher<true,TID>::transferForward(TComm<DIMENSION,TTypes>& comm, typename TTypes<TID>::FwdType& rData)
   {
   }

   template <Dimensions::Transform::Id TID> class EndDispatcher<false,TID>
   {
      public:
         /**
          * @brief Receive forward the data
          */
         template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static typename TTypes<TID>::FwdType& receiveForward(TComm<DIMENSION,TTypes>& comm);

         /**
          * @brief transfer forward data
          */
         template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> static void transferForward(TComm<DIMENSION,TTypes>& comm, typename TTypes<TID>::FwdType& rData);

   };

   template <Dimensions::Transform::Id TID> template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> typename TTypes<TID>::FwdType& EndDispatcher<false,TID>::receiveForward(TComm<DIMENSION,TTypes>& comm)
   {
      typename TTypes<TID>::FwdType &rData = comm.template converter<Dimensions::Transform::jump<TID,1>::id>().getFwd(comm.template storage<TID>());

      return rData;
   }

   template <Dimensions::Transform::Id TID> template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes, template <Dimensions::Type, template<Dimensions::Transform::Id> class > class TComm> void EndDispatcher<false,TID>::transferForward(TComm<DIMENSION,TTypes>& comm, typename TTypes<TID>::FwdType& rData)
   {
    // Convert data
    comm.template converter<Dimensions::Transform::jump<TID,1>::id>().convertFwd(rData, comm.template storage<Dimensions::Transform::jump<TID,1>::id>());

    // Free input data
    comm.template storage<TID>().freeFwd(rData);
   }
}
}

#endif // ENDDISPATCHER_HPP
