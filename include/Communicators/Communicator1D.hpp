/** \file Communicator1D.hpp
 *  \brief Implementation of a 1D communicator
 */

#ifndef COMMUNICATOR1D_HPP
#define COMMUNICATOR1D_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "StorageProviders/StoragePairProvider.hpp"
#include "Communicators/CommunicatorBase.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of 1D communicator
    *
    * @tparam TFwd1D Data type for the forward 1D transform
    * @tparam TBwd1D Data type for the backward 1D transform
    */ 
   template <typename TFwd1D, typename TBwd1D> class Communicator1D: public CommunicatorBase
   {
      public:
         /// Typedef for the forward data type in 1D
         typedef TFwd1D Fwd1DType;

         /// Typdeef for the backward data type in 1D
         typedef TBwd1D Bwd1DType;

         /**
         * @brief Very basic constructor
         */
         Communicator1D();

         /**
         * @brief Destructor
         */
         virtual ~Communicator1D();

         /**
          * @brief Initialise the coordinator
          *
          * @param setupFwd1D Setup object for the forward 1D type
          * @param setupBwd1D Setup object for the backward 1D type
          */
         void init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D);

         /**
          * @brief Get the storage provider for the first transform
          */
         StoragePairProvider<TFwd1D, TBwd1D>&  storage1D();

         /**
          * @brief Recover the TFwd1D from storage
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         TFwd1D&  receiveFwd1D();

         /**
          * @brief Transfer TFwd1D to next step
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         void transferFwd1D(TFwd1D& rData);

         /**
          * @brief Hold starting spectral TBwd1D
          */
         void holdSpectral(TBwd1D& rData);

         /**
          * @brief Free the starting spectral TBwd1D
          *
          * This routine does nothing but is used to make sure storage doesn't end up in storage provider loop
          */
         void freeSpectral(TBwd1D& rData);

         /**
          * @brief Provide physical storage
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         TFwd1D&  providePhysical();

         /**
          * @brief Hold final physical TFwd1D
          */
         void holdPhysical(TFwd1D& rData);
         
      protected:

      private:
         /**
          * @brief The storage provider for the first dimension
          */
         StoragePairProvider<TFwd1D, TBwd1D>  mStorage1D;
   };

   template <typename TFwd1D, typename TBwd1D> inline StoragePairProvider<TFwd1D, TBwd1D>&  Communicator1D<TFwd1D, TBwd1D>::storage1D()
   {
      return this->mStorage1D;
   }

   template <typename TFwd1D, typename TBwd1D> Communicator1D<TFwd1D, TBwd1D>::Communicator1D()
   {
   }

   template <typename TFwd1D, typename TBwd1D> Communicator1D<TFwd1D, TBwd1D>::~Communicator1D()
   {
   }

   template <typename TFwd1D, typename TBwd1D> void Communicator1D<TFwd1D, TBwd1D>::init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D)
   {
      // Initialise first dimension storage
      this->mStorage1D.init(setupFwd1D, setupBwd1D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem1D = this->mStorage1D.requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem1D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM1D, mem1D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TFwd1D, typename TBwd1D> TFwd1D& Communicator1D<TFwd1D, TBwd1D>::receiveFwd1D()
   {
      return this->storage1D().recoverFwd();
   }

   template <typename TFwd1D, typename TBwd1D> void Communicator1D<TFwd1D, TBwd1D>::transferFwd1D(TFwd1D& rData)
   {
   }

   template <typename TFwd1D, typename TBwd1D> void Communicator1D<TFwd1D, TBwd1D>::holdSpectral(TBwd1D& rData)
   {
      // Hold the input data
      this->storage1D().holdBwd(rData);
   }

   template <typename TFwd1D, typename TBwd1D> void Communicator1D<TFwd1D, TBwd1D>::freeSpectral(TBwd1D& rData)
   {
      // 
      // DO NOTHING BUT MAKE SURE IT IS THE CASE
      //
   }

   template <typename TFwd1D, typename TBwd1D> TFwd1D& Communicator1D<TFwd1D, TBwd1D>::providePhysical()
   {
      return this->storage1D().provideFwd();
   }

   template <typename TFwd1D, typename TBwd1D> void Communicator1D<TFwd1D, TBwd1D>::holdPhysical(TFwd1D& rData)
   {
      // Hold the input data
      this->storage1D().holdFwd(rData);
   }

}
}

#endif // COMMUNICATOR1D_HPP
