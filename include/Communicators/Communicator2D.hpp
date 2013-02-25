/** \file Communicator2D.hpp
 *  \brief Implementation of a 2D communicator
 */

#ifndef COMMUNICATOR2D_HPP
#define COMMUNICATOR2D_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "TypeSelectors/IndexConverterSelector.hpp"
#include "StorageProviders/StoragePairProvider.hpp"
#include "Communicators/Communicator1D.hpp"
#include "Communicators/Converters/IConverter.hpp"
#include "Communicators/Converters/SerialConverter.hpp"
#ifdef GEOMHDISCC_MPI
   #include "Communicators/Converters/MpiConverter.hpp"
#endif // GEOMHDISCC_MPI

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of 2D communicator
    *
    * @tparam TFwd1D Data type for the forward 1D transform
    * @tparam TBwd1D Data type for the backward 1D transform
    * @tparam TFwd2D Data type for the forward 2D transform
    * @tparam TBwd2D Data type for the backward 2D transform
    */ 
   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> class Communicator2D: public Communicator1D<TFwd1D, TBwd1D>
   {
      public:
         /// Typedef for the forward data type in 1D
         typedef TFwd1D Fwd1DType;

         /// Typedef for the backward data type in 1D
         typedef TBwd1D Bwd1DType;

         /// Typedef for the forward data type in 2D
         typedef TFwd2D Fwd2DType;

         /// Typedef for the backward data type in 2D
         typedef TBwd2D Bwd2DType;

         /**
         * @brief Very basic constructor
         */
         Communicator2D();

         /**
         * @brief Destructor
         */
         virtual ~Communicator2D();

         /**
          * @brief Initialise the coordinator
          *
          * @param setupFwd1D Setup object for the forward 1D type
          * @param setupBwd1D Setup object for the backward 1D type
          * @param setupFwd2D Setup object for the forward 2D type
          * @param setupBwd2D Setup object for the backward 2D type
          */
         void init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D, const typename TFwd2D::SetupType& setupFwd2D, const typename TBwd2D::SetupType& setupBwd2D);

         /**
          * @brief Initialise the converter
          *
          * @param spRes      Shared resolution information
          * @param packs1DFwd Packs information for first forward exchange
          * @param packs1DBwd Packs information for first backward exchange
          * @param split      Location where the MPI splitting takes place
          */
         void initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, Splitting::Locations::Id split);

         /**
          * @brief Get the storage provider for the second transform
          */
         StoragePairProvider<TFwd2D, TBwd2D>&  storage2D();

         /**
          * @brief Receive the TFwd1D transfered as an TBwd2D
          */
         TFwd1D&  receiveFwd1D();

         /**
          * @brief Recover the TFwd2D from storage
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         TFwd2D&  receiveFwd2D();

         /**
          * @brief Receive the TBwd2D transfered as an TFwd1D
          */
         TBwd2D&  receiveBwd2D();

         /**
          * @brief Get the first to second transform converter
          */
         IConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D>&  converter2D();

         /**
          * @brief Transfer TFwd1D to next step
          */
         void transferFwd1D(TFwd1D& rData);

         /**
          * @brief Transfer TFwd2D to next step
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         void transferFwd2D(TFwd2D& rData);

         /**
          * @brief Transfer TBwd2D to next step
          */
         void transferBwd2D(TBwd2D& rData);

         /**
          * @brief Provide physical storage
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         TFwd2D&  providePhysical();

         /**
          * @brief Hold final physical TFwd2D
          */
         void holdPhysical(TFwd2D& rData);
         
      protected:
         /**
          * @brief Converter between first and second transform
          */
         SharedPtrMacro<IConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D> > mspConverter2D;

      private:
         /**
          * @brief The storage provider for the second dimension
          */
         StoragePairProvider<TFwd2D, TBwd2D>  mStorage2D;
   };

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> inline StoragePairProvider<TFwd2D, TBwd2D>&  Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::storage2D()
   {
      return this->mStorage2D;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> inline IConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D>& Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::converter2D()
   {
      return *this->mspConverter2D;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::Communicator2D()
      : Communicator1D<TFwd1D, TBwd1D>()
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::~Communicator2D()
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D, const typename TFwd2D::SetupType& setupFwd2D, const typename TBwd2D::SetupType& setupBwd2D)
   {
      // Initialise 1D storage
      Communicator1D<TFwd1D,TBwd1D>::init(setupFwd1D, setupBwd1D);

      // Initialise second dimension storage
      this->mStorage2D.init(setupFwd2D, setupBwd2D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem2D = this->mStorage2D.requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem2D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM2D, mem2D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, Splitting::Locations::Id split)
   {
      #ifdef GEOMHDISCC_MPI
         if(split == Splitting::Locations::FIRST)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<MpiConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D> > spConv(new MpiConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D>());
            spConv->init(spRes, 0, this->storage1D().rFTmps(), this->storage2D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(2);

            // Set communicattion buffers' pointers
            this->mspConverter2D->setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));

            // Get maximum number of packs
            int max2D = std::max(packs1DFwd(packs1DFwd.size()-1), packs1DBwd(packs1DBwd.size()-1));

            // Allocate first 2D buffers
            this->allocateBuffers(0, this->mspConverter2D->fwdSizes(), max2D);

            // Allocate second 2D buffers
            this->allocateBuffers(1, this->mspConverter2D->bwdSizes(), max2D);
         } else if(split == Splitting::Locations::SECOND)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<SerialConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type> > spConv(new SerialConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>());

            // initialise the converter
            spConv->init(spRes, 0);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(2);
         } else if(split == Splitting::Locations::BOTH)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<MpiConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D> > spConv(new MpiConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D>());
            spConv->init(spRes, 0, this->storage1D().rFTmps(), this->storage2D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(3);

            // Set communicattion buffers' pointers
            this->mspConverter2D->setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));
         }
      #else 
         // Initialise 2D serial converter
         SharedPtrMacro<SerialConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type> > spConv(new SerialConverter<TFwd1D, TBwd1D, TFwd2D, TBwd2D, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>());

         // initialise the converter
         spConv->init(spRes, 0);

         this->mspConverter2D = spConv;
      #endif // GEOMHDISCC_MPI

      // If only one dimension is split
      if(split != Splitting::Locations::BOTH)
      {
         // Setup converter
         this->mspConverter2D->setup();
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         // Do (MPI) storage profiling on 2D converter
         this->mspConverter2D->profileStorage();
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> TFwd1D&  Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::receiveFwd1D()
   {
      TFwd1D &rData = this->mspConverter2D->getFwd(this->storage1D());

      return rData;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> TFwd2D&  Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::receiveFwd2D()
   {
      return this->storage2D().recoverFwd();
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> TBwd2D&  Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::receiveBwd2D()
   {
      TBwd2D &rData = this->mspConverter2D->getBwd(this->storage2D());

      return rData;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::transferFwd1D(TFwd1D& rData)
   {
      // Convert data
      this->mspConverter2D->convertFwd(rData, this->storage2D());

      // Free input data
      this->storage1D().freeFwd(rData);
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::transferFwd2D(TFwd2D& rData)
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::transferBwd2D(TBwd2D& rData)
   {
      // Convert data
      this->mspConverter2D->convertBwd(rData, this->storage1D());

      // Free input data
      this->storage2D().freeBwd(rData);
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> TFwd2D&  Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::providePhysical()
   {
      return this->storage2D().provideFwd();
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D> void Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>::holdPhysical(TFwd2D& rData)
   {
      // Hold the input data
      this->storage2D().holdFwd(rData);
   }

}
}

#endif // COMMUNICATOR2D_HPP
