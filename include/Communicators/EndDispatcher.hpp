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
         template <typename TData, typename TComm> static TData& receiveForward(TComm& comm);

         /**
          * @brief transfer forward data
          */
         template <typename TComm, typename TData> static void transferForward(TComm& comm, TData& rData);

   };

   template <Dimensions::Transform::Id TID> template <typename TData, typename TComm> TData& EndDispatcher<true,TID>::receiveForward(TComm& comm)
   {
      return comm.template storage<TID>().recoverFwd();
   }

   template <Dimensions::Transform::Id TID> template <typename TComm, typename TData> void EndDispatcher<true,TID>::transferForward(TComm& comm, TData& rData)
   {
   }

   template <Dimensions::Transform::Id TID> class EndDispatcher<false,TID>
   {
      public:
         /**
          * @brief Receive forward the data
          */
         template <typename TData, typename TComm> static TData& receiveForward(TComm& comm);

         /**
          * @brief transfer forward data
          */
         template <typename TComm, typename TData> static void transferForward(TComm& comm, TData& rData);

   };

   template <Dimensions::Transform::Id TID> template <typename TData, typename TComm> TData& EndDispatcher<false,TID>::receiveForward(TComm& comm)
   {
      TData &rData = comm.template converter<Dimensions::Transform::jump<TID,1>::id>().getFwd(comm.template storage<TID>());

      return rData;
   }

   template <Dimensions::Transform::Id TID> template <typename TComm, typename TData> void EndDispatcher<false,TID>::transferForward(TComm& comm, TData& rData)
   {
    // Convert data
    comm.template converter<Dimensions::Transform::jump<TID,1>::id>().convertFwd(rData, comm.template storage<Dimensions::Transform::jump<TID,1>::id>());

    // Free input data
    comm.template storage<TID>().freeFwd(rData);
   }
}
}

#endif // ENDDISPATCHER_HPP
