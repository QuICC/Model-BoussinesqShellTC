/** \file Communicator.hpp
 *  \brief Implementation of a 3D communicator
 *
 *  \mhdBug Needs test
 */

#ifndef COMMUNICATOR_HPP
#define COMMUNICATOR_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/Splitting.hpp"
#include "Resolutions/Resolution.hpp"
#include "Communicators/CommunicatorStorage.hpp"
#include "Communicators/EndDispatcher.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of 3D communicator
    */ 
   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> class Communicator: public CommunicatorStorage<DIMENSION,TTypes>
   {
      public:
         /**
         * @brief Very basic constructor
         */
         Communicator();

         /**
         * @brief Destructor
         */
         ~Communicator();

         /**
          * @brief Initialise the coordinator
          *
          * @param setupFwd1D Setup object for the forward 1D type
          * @param setupBwd1D Setup object for the backward 1D type
          */
         void init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D);

         /**
          * @brief Initialise the coordinator
          *
          * @param setupFwd1D Setup object for the forward 1D type
          * @param setupBwd1D Setup object for the backward 1D type
          * @param setupFwd2D Setup object for the forward 2D type
          * @param setupBwd2D Setup object for the backward 2D type
          */
         void init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D, const typename TTypes<Dimensions::Transform::TRA2D>::FwdType::SetupType& setupFwd2D, const typename TTypes<Dimensions::Transform::TRA2D>::BwdType::SetupType& setupBwd2D);

         /**
          * @brief Initialise the coordinator
          *
          * @param setupFwd1D Setup object for the forward 1D type
          * @param setupBwd1D Setup object for the backward 1D type
          * @param setupFwd2D Setup object for the forward 2D type
          * @param setupBwd2D Setup object for the backward 2D type
          * @param setupFwd3D Setup object for the forward 3D type
          * @param setupBwd3D Setup object for the backward 3D type
          */
         void init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D, const typename TTypes<Dimensions::Transform::TRA2D>::FwdType::SetupType& setupFwd2D, const typename TTypes<Dimensions::Transform::TRA2D>::BwdType::SetupType& setupBwd2D, const typename TTypes<Dimensions::Transform::TRA3D>::FwdType::SetupType& setupFwd3D, const typename TTypes<Dimensions::Transform::TRA3D>::BwdType::SetupType& setupBwd3D);

         /**
          * @brief Initialise the 1D/2D converter
          *
          * @param spRes      Shared resolution information
          * @param packs1DFwd Packs information for first forward exchange
          * @param packs1DBwd Packs information for first backward exchange
          * @param split      Location where the MPI splitting takes place
          */
         void initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, Splitting::Locations::Id split);

         /**
          * @brief Initialise the 2D/3D converter
          *
          * @param spRes      Shared resolution information
          * @param packs1DFwd Packs information for first forward exchange
          * @param packs1DBwd Packs information for first backward exchange
          * @param packs2DFwd Packs information for second forward exchange
          * @param packs2DBwd Packs information for second backward exchange
          * @param split      Location where the MPI splitting takes place
          */
         void initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, const ArrayI& packs2DFwd, const ArrayI& packs2DBwd, Splitting::Locations::Id split);

         /**
          * @brief Transfer forward data to next step
          */
         template <Dimensions::Transform::Id TID> void transferForward(typename TTypes<TID>::FwdType& rData);

         /**
          * @brief Transfer backward data to next step
          */
         template <Dimensions::Transform::Id TID> void transferBackward(typename TTypes<TID>::BwdType& rData);

         /**
          * @brief Receive forward the data
          */
         template <Dimensions::Transform::Id TID> typename TTypes<TID>::FwdType& receiveForward();
         /**
          * @brief Receive backward the data
          */
         template <Dimensions::Transform::Id TID> typename TTypes<TID>::BwdType& receiveBackward();

         /**
          * @brief Provide physical storage
          */
         typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>::FwdType&  providePhysical();

         /**
          * @brief Hold final physical storage
          */
         void holdPhysical(typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>::FwdType& rData);

         /**
          * @brief Dealias the spectral data
          */
         void dealiasSpectral(const typename TTypes<Dimensions::Transform::TRA1D>::BwdType& rData);
         
      protected:

      private:
         template <bool COND, Dimensions::Transform::Id TID> friend class EndDispatcher;

   };

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> Communicator<DIMENSION,TTypes>::Communicator()
   {
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> Communicator<DIMENSION,TTypes>::~Communicator()
   {
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> typename TTypes<TID>::FwdType&  Communicator<DIMENSION,TTypes>::receiveForward()
   {
      Debug::StaticAssert< static_cast<int>(TID) <= static_cast<int>(DIMENSION) >();

      return EndDispatcher<(static_cast<int>(DIMENSION) == static_cast<int>(TID)), TID>::receiveForward(*this);
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> typename TTypes<TID>::BwdType&  Communicator<DIMENSION,TTypes>::receiveBackward()
   {
      Debug::StaticAssert< (static_cast<int>(TID) > 0) && (static_cast<int>(TID) <= static_cast<int>(DIMENSION)) >();

      typename TTypes<TID>::BwdType &rData = this->template converter<TID>().getBwd(this->template storage<TID>());

      return rData;
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> void Communicator<DIMENSION,TTypes>::transferForward(typename TTypes<TID>::FwdType& rData)
   {
      Debug::StaticAssert< (static_cast<int>(TID) <= static_cast<int>(DIMENSION)) >();

      EndDispatcher<(static_cast<int>(DIMENSION) == static_cast<int>(TID)), TID>::transferForward(*this, rData);
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> template <Dimensions::Transform::Id TID> void Communicator<DIMENSION,TTypes>::transferBackward(typename TTypes<TID>::BwdType& rData)
   {
      Debug::StaticAssert< (static_cast<int>(TID) > 0) && (static_cast<int>(TID) <= static_cast<int>(DIMENSION)) >();

      // Convert data
      this->template converter<TID>().convertBwd(rData, this->template storage<Dimensions::Transform::jump<TID,-1>::id>());

      // Free input data
      this->template storage<TID>().freeBwd(rData);
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>::FwdType&  Communicator<DIMENSION,TTypes>::providePhysical()
   {
      return this->template storage<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>().provideFwd();
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::holdPhysical(typename TTypes<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>::FwdType& rData)
   {
      // Hold the input data
      this->template storage<static_cast<Dimensions::Transform::Id>(static_cast<int>(DIMENSION))>().holdFwd(rData);
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::dealiasSpectral(const typename TTypes<Dimensions::Transform::TRA1D>::BwdType& rInData)
   {
      // Get dealiased storage scalar
      typename TTypes<Dimensions::Transform::TRA1D>::BwdType &rOutData = this->template storage<Dimensions::Transform::TRA1D>().provideBwd();

      // Dealias the data
      rOutData.rData().topRows(rInData.data().rows()) = rInData.data(); 

      // Hold the input data
      this->template storage<Dimensions::Transform::TRA1D>().holdBwd(rOutData);
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D)
   {
      Debug::StaticAssert< static_cast<int>(DIMENSION) >= 0 >();

      // Initialise storage for 1D
      this->template storage<Dimensions::Transform::TRA1D>().init(setupFwd1D, setupBwd1D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem1D = this->template storage<Dimensions::Transform::TRA1D>().requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem1D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM3D, mem1D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D, const typename TTypes<Dimensions::Transform::TRA2D>::FwdType::SetupType& setupFwd2D, const typename TTypes<Dimensions::Transform::TRA2D>::BwdType::SetupType& setupBwd2D)
   {
      Debug::StaticAssert< (static_cast<int>(DIMENSION) > 0) >();

      // Initialise 1D storage
      this->init(setupFwd1D, setupBwd1D);

      // Initialise storage for 2D
      this->template storage<Dimensions::Transform::TRA2D>().init(setupFwd2D, setupBwd2D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem2D = this->template storage<Dimensions::Transform::TRA2D>().requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem2D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM3D, mem2D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::init(const typename TTypes<Dimensions::Transform::TRA1D>::FwdType::SetupType& setupFwd1D, const typename TTypes<Dimensions::Transform::TRA1D>::BwdType::SetupType& setupBwd1D, const typename TTypes<Dimensions::Transform::TRA2D>::FwdType::SetupType& setupFwd2D, const typename TTypes<Dimensions::Transform::TRA2D>::BwdType::SetupType& setupBwd2D, const typename TTypes<Dimensions::Transform::TRA3D>::FwdType::SetupType& setupFwd3D, const typename TTypes<Dimensions::Transform::TRA3D>::BwdType::SetupType& setupBwd3D)
   {
      Debug::StaticAssert< (static_cast<int>(DIMENSION) > 1) >();

      // Initialise 2D storage
      this->init(setupFwd1D, setupBwd1D, setupFwd2D, setupBwd2D);

      // Initialise storage for 3D
      this->template storage<Dimensions::Transform::TRA3D>().init(setupFwd3D, setupBwd3D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem3D = this->template storage<Dimensions::Transform::TRA3D>().requiredStorage();
         StorageProfilerMacro_update(Debug::StorageProfiler::TEMPORARIES, mem3D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(Debug::StorageProfiler::TRANSFORM3D, mem3D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, Splitting::Locations::Id split)
   {
      #ifdef GEOMHDISCC_MPI
         if(split == Splitting::Locations::FIRST)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<MpiConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType> > spConv(new MpiConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType>());
            spConv->init(spRes, 0, this->storage1D().rFTmps(), this->storage2D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(2);

            // Set communicattion buffers' pointers
            this->template converter<Dimensions::Transform::TRA2D>().setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));

            // Get maximum number of packs
            int max2D = std::max(packs1DFwd(packs1DFwd.size()-1), packs1DBwd(packs1DBwd.size()-1));

            // Allocate first 2D buffers
            this->allocateBuffers(0, this->template converter<Dimensions::Transform::TRA2D>()->fwdSizes(), max2D);

            // Allocate second 2D buffers
            this->allocateBuffers(1, this->template converter<Dimensions::Transform::TRA2D>()->bwdSizes(), max2D);
         } else if(split == Splitting::Locations::SECOND)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<SerialConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type> > spConv(new SerialConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>());

            // initialise the converter
            spConv->init(spRes, 0);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(2);
         } else if(split == Splitting::Locations::BOTH)
         {
            // Initialise 2D serial converter
            SharedPtrMacro<MpiConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType> > spConv(new MpiConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType>());
            spConv->init(spRes, 0, this->storage1D().rFTmps(), this->storage2D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter2D = spConv;

            // Create the communication buffers
            this->createBuffers(3);

            // Set communicattion buffers' pointers
            this->template converter<Dimensions::Transform::TRA2D>()->setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));
         }
      #else 
         // Initialise 2D serial converter
         SharedPtrMacro<SerialConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type> > spConv(new SerialConverter<typename TTypes<Dimensions::Transform::TRA1D>::FwdType, typename TTypes<Dimensions::Transform::TRA1D>::BwdType, typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA2D>::Type>());

         // initialise the converter
         spConv->init(spRes, 0);

         this->mspConverter2D = spConv;
      #endif // GEOMHDISCC_MPI

      // If only one dimension is split
      if(split != Splitting::Locations::BOTH)
      {
         // Setup converter
         this->template converter<Dimensions::Transform::TRA2D>().setup();
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         // Do (MPI) storage profiling on 2D converter
         this->template converter<Dimensions::Transform::TRA2D>().profileStorage();
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <Dimensions::Type DIMENSION, template <Dimensions::Transform::Id> class TTypes> void Communicator<DIMENSION,TTypes>::initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, const ArrayI& packs2DFwd, const ArrayI& packs2DBwd, Splitting::Locations::Id split)
   {
      // Initialise 1D/2D converter
      this->initConverter(spRes, packs1DFwd, packs1DBwd, split);
     
      #ifdef GEOMHDISCC_MPI
         if(split == Splitting::Locations::FIRST)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<SerialConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA3D>::Type> > spConv(new SerialConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA3D>::Type>());

            // Initialise the converter
            spConv->init(spRes, 1);

            this->mspConverter3D = spConv;

         } else if(split == Splitting::Locations::SECOND)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<MpiConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType> > spConv(new MpiConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType>());
            spConv->init(spRes, 1, this->storage2D().rFTmps(), this->storage3D().rBTmps(), packs2DFwd, packs2DBwd);

            this->mspConverter3D = spConv;

            // Set communicattion buffers' pointers
            this->template converter<Dimensions::Transform::TRA3D>().setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));

            // Get maximum number of packs
            int max3D = std::max(packs2DFwd(packs2DFwd.size()-1), packs2DBwd(packs2DBwd.size()-1));

            // Allocate first 3D buffers
            this->allocateBuffers(0, this->template converter<Dimensions::Transform::TRA3D>().fwdSizes(), max3D);

            // Allocate second 3D buffers
            this->allocateBuffers(1, this->template converter<Dimensions::Transform::TRA3D>().bwdSizes(), max3D);

         } else if(split == Splitting::Locations::BOTH)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<MpiConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType> > spConv(new MpiConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType>());
            spConv->init(spRes, 1, this->storage2D().rFTmps(), this->storage3D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter3D = spConv;

            // Set communicattion buffers' pointers
            this->template converter<Dimensions::Transform::TRA3D>().setBuffers(this->mBuffers.at(2), this->mBuffers.at(0));

            // Get maximum number of packs
            int max2D = std::max(packs1DFwd(packs1DFwd.size()-1), packs1DBwd(packs1DBwd.size()-1));

            // Get maximum number of packs
            int max3D = std::max(packs2DFwd(packs2DFwd.size()-1), packs2DBwd(packs2DBwd.size()-1));

            // Allocate shared buffers
            this->allocateBuffers(0, this->template converter<Dimensions::Transform::TRA2D>().fwdSizes(), max2D, this->template converter<Dimensions::Transform::TRA3D>().bwdSizes(), max3D);

            // Allocate 2D buffers
            this->allocateBuffers(1, this->template converter<Dimensions::Transform::TRA2D>().wdSizes(), max2D);

            // Allocate 3D buffers
            this->allocateBuffers(2, this->template converter<Dimensions::Transform::TRA3D>().fwdSizes(), max3D);
         }
      #else 
         // Initialise serial 3D converter
         SharedPtrMacro<SerialConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA3D>::Type> > spConv(new SerialConverter<typename TTypes<Dimensions::Transform::TRA2D>::FwdType, typename TTypes<Dimensions::Transform::TRA2D>::BwdType, typename TTypes<Dimensions::Transform::TRA3D>::FwdType, typename TTypes<Dimensions::Transform::TRA3D>::BwdType, IndexConverterSelector<Dimensions::Transform::TRA3D>::Type>());

         // Initialise the converter
         spConv->init(spRes, 1);

         this->mspConverter3D = spConv;
      #endif // GEOMHDISCC_MPI

      // Setup converter
      this->template converter<Dimensions::Transform::TRA3D>().setup();

      // If both dimensions are split
      if(split == Splitting::Locations::BOTH)
      {
         this->template converter<Dimensions::Transform::TRA2D>().setup();
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         // Do (MPI) storage profiling on 3D converter
         this->template converter<Dimensions::Transform::TRA3D>().profileStorage();
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

}
}

#endif // COMMUNICATOR_HPP
