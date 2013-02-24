/** \file Communicator3D.hpp
 *  \brief Implementation of a 3D communicator
 */

#ifndef COMMUNICATOR3D_HPP
#define COMMUNICATOR3D_HPP

// Configuration includes
//
#include "Base/PrepMacros/SmartPointerMacro.h"
#include "Base/PrepMacros/IndexConverterMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Enums/Splittings.hpp"
#include "Base/Resolutions/Resolution.hpp"
#include "Base/StorageProviders/StoragePairProvider.hpp"
#include "Base/Communicators/Communicator2D.hpp"
#include "Base/Communicators/Converters/ConverterBase.hpp"
#include "Base/Communicators/Converters/SerialConverter.hpp"
#ifdef GEOMHDISCC_MPI
   #include "Base/Communicators/Converters/MpiConverter.hpp"
#endif // GEOMHDISCC_MPI

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of 3D communicator
    *
    * @tparam TFwd1D Data type for the forward 1D transform
    * @tparam TBwd1D Data type for the backward 1D transform
    * @tparam TFwd2D Data type for the forward 2D transform
    * @tparam TBwd2D Data type for the backward 2D transform
    * @tparam TFwd3D Data type for the forward 3D transform
    * @tparam TBwd3D Data type for the backward 3D transform
    */ 
   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> class Communicator3D: public Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>
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

         /// Typedef for the forward data type in 3D
         typedef TFwd3D Fwd3DType;

         /// Typedef for the backward data type in 3D
         typedef TBwd3D Bwd3DType;

         /**
         * @brief Very basic constructor
         */
         Communicator3D();

         /**
         * @brief Destructor
         */
         virtual ~Communicator3D();

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
         void init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D, const typename TFwd2D::SetupType& setupFwd2D, const typename TBwd2D::SetupType& setupBwd2D, const typename TFwd3D::SetupType& setupFwd3D, const typename TBwd3D::SetupType& setupBwd3D);

         /**
          * @brief Initialise the converter
          *
          * @param spRes      Shared resolution information
          * @param packs1DFwd Packs information for first forward exchange
          * @param packs1DBwd Packs information for first backward exchange
          * @param packs2DFwd Packs information for second forward exchange
          * @param packs2DBwd Packs information for second backward exchange
          * @param split      Location where the MPI splitting takes place
          */
         void initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, const ArrayI& packs2DFwd, const ArrayI& packs2DBwd, Splittings::Locations::Id split);

         /**
          * @brief Get the storage provider for the third transform
          */
         StoragePairProvider<TFwd3D, TBwd3D>&  storage3D();

         /**
          * @brief Receive the TFwd2D transfered as an TBwd3D
          */
         TFwd2D&  receiveFwd2D();

         /**
          * @brief Recover the TFwd3D from storage
          *
          * This routine is required to generalise to 1D, 2D, 3D. 
          * It has to be overloaded in higher dimensions.
          */
         TFwd3D&  receiveFwd3D();

         /**
          * @brief Receive the TBwd3D transfered as an TFwd2D
          */
         TBwd3D&  receiveBwd3D();

         /**
          * @brief Get the second to third transform converter
          */
         ConverterBase<TFwd2D, TBwd2D, TFwd3D, TBwd3D>&  converter3D();

         /**
          * @brief Transfer TFwd2D to next step
          */
         void transferFwd2D(TFwd2D& rData);

         /**
          * @brief Transfer TFwd3D to next step
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         void transferFwd3D(TFwd3D& rData);

         /**
          * @brief Transfer TBwd3D to next step
          */
         void transferBwd3D(TBwd3D& rData);

         /**
          * @brief Provide physical storage
          *
          * This routine is required to generalise to 1D, 2D, 3D.
          * It has to be overloaded in higher dimensions.
          */
         TFwd3D&  providePhysical();

         /**
          * @brief Hold final physical TFwd3D
          */
         void holdPhysical(TFwd3D& rData);
         
      protected:
         /**
          * @brief Converter between second and third transform
          */
         SharedPtrMacro<ConverterBase<TFwd2D, TBwd2D, TFwd3D, TBwd3D> > mspConverter3D;

      private:
         /**
          * @brief The storage provider for the third dimension
          */
         StoragePairProvider<TFwd3D, TBwd3D>  mStorage3D;
   };

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> inline StoragePairProvider<TFwd3D, TBwd3D>&  Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::storage3D()
   {
      return this->mStorage3D;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> inline ConverterBase<TFwd2D, TBwd2D, TFwd3D, TBwd3D>& Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::converter3D()
   {
      return *this->mspConverter3D;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::Communicator3D()
      : Communicator2D<TFwd1D, TBwd1D, TFwd2D, TBwd2D>()
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::~Communicator3D()
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::init(const typename TFwd1D::SetupType& setupFwd1D, const typename TBwd1D::SetupType& setupBwd1D, const typename TFwd2D::SetupType& setupFwd2D, const typename TBwd2D::SetupType& setupBwd2D, const typename TFwd3D::SetupType& setupFwd3D, const typename TBwd3D::SetupType& setupBwd3D)
   {
      // Initialise 2D storage
      Communicator2D<TFwd1D,TBwd1D, TFwd2D,TBwd2D>::init(setupFwd1D, setupBwd1D, setupFwd2D, setupBwd2D);

      // Initialise third dimension storage
      this->mStorage3D.init(setupFwd3D, setupBwd3D);

      #ifdef GEOMHDISCC_STORAGEPROFILE
         MHDFloat mem3D = this->mStorage3D.requiredStorage();
         StorageProfilerMacro_update(StorageProfiler::TEMPORARIES, mem3D);

         #ifdef GEOMHDISCC_STORAGEPROFILER_DETAILED
            StorageProfilerMacro_update(StorageProfiler::TRANSFORM3D, mem3D);
         #endif // GEOMHDISCC_STORAGEPROFILER_DETAILED
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::initConverter(SharedResolution spRes, const ArrayI& packs1DFwd, const ArrayI& packs1DBwd, const ArrayI& packs2DFwd, const ArrayI& packs2DBwd, Splittings::Locations::Id split)
   {
      // Initialise 2D converter
      Communicator2D<TFwd1D,TBwd1D, TFwd2D,TBwd2D>::initConverter(spRes, packs1DFwd, packs1DBwd, split);
     
      #ifdef GEOMHDISCC_MPI
         if(split == Splittings::Locations::FIRST)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<SerialConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D, IndexConverter3DMacro> > spConv(new SerialConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D, IndexConverter3DMacro>());

            // Initialise the converter
            spConv->init(spRes, 1);

            this->mspConverter3D = spConv;

         } else if(split == Splittings::Locations::SECOND)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<MpiConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D> > spConv(new MpiConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D>());
            spConv->init(spRes, 1, this->storage2D().rFTmps(), this->storage3D().rBTmps(), packs2DFwd, packs2DBwd);

            this->mspConverter3D = spConv;

            // Set communicattion buffers' pointers
            this->mspConverter3D->setBuffers(this->mBuffers.at(0), this->mBuffers.at(1));

            // Get maximum number of packs
            int max3D = std::max(packs2DFwd(packs2DFwd.size()-1), packs2DBwd(packs2DBwd.size()-1));

            // Allocate first 3D buffers
            this->allocateBuffers(0, this->mspConverter3D->fwdSizes(), max3D);

            // Allocate second 3D buffers
            this->allocateBuffers(1, this->mspConverter3D->bwdSizes(), max3D);

         } else if(split == Splittings::Locations::BOTH)
         {
            // Initialise serial 3D converter
            SharedPtrMacro<MpiConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D> > spConv(new MpiConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D>());
            spConv->init(spRes, 1, this->storage2D().rFTmps(), this->storage3D().rBTmps(), packs1DFwd, packs1DBwd);

            this->mspConverter3D = spConv;

            // Set communicattion buffers' pointers
            this->mspConverter3D->setBuffers(this->mBuffers.at(2), this->mBuffers.at(0));

            // Get maximum number of packs
            int max2D = std::max(packs1DFwd(packs1DFwd.size()-1), packs1DBwd(packs1DBwd.size()-1));

            // Get maximum number of packs
            int max3D = std::max(packs2DFwd(packs2DFwd.size()-1), packs2DBwd(packs2DBwd.size()-1));

            // Allocate shared buffers
            this->allocateBuffers(0, this->mspConverter2D->fwdSizes(), max2D, this->mspConverter3D->bwdSizes(), max3D);

            // Allocate 2D buffers
            this->allocateBuffers(1, this->mspConverter2D->bwdSizes(), max2D);

            // Allocate 3D buffers
            this->allocateBuffers(2, this->mspConverter3D->fwdSizes(), max3D);
         }
      #else 
         // Initialise serial 3D converter
         SharedPtrMacro<SerialConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D, IndexConverter3DMacro> > spConv(new SerialConverter<TFwd2D, TBwd2D, TFwd3D, TBwd3D, IndexConverter3DMacro>());

         // Initialise the converter
         spConv->init(spRes, 1);

         this->mspConverter3D = spConv;
      #endif // GEOMHDISCC_MPI

      // Setup converter
      this->mspConverter3D->setup();

      // If both dimensions are split
      if(split == Splittings::Locations::BOTH)
      {
         this->mspConverter2D->setup();
      }

      #ifdef GEOMHDISCC_STORAGEPROFILE
         // Do (MPI) storage profiling on 3D converter
         this->mspConverter3D->profileStorage();
      #endif // GEOMHDISCC_STORAGEPROFILE
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> TFwd2D&  Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::receiveFwd2D()
   {
      TFwd2D &rData = this->mspConverter3D->getFwd(this->storage2D());

      return rData;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> TFwd3D&  Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::receiveFwd3D()
   {
      return this->storage3D().recoverFwd();
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> TBwd3D&  Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::receiveBwd3D()
   {
      TBwd3D &rData = this->mspConverter3D->getBwd(this->storage3D());

      return rData;
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::transferFwd2D(TFwd2D& rData)
   {
      // Convert data
      this->mspConverter3D->convertFwd(rData, this->storage3D());

      // Free input data
      this->storage2D().freeFwd(rData);
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::transferFwd3D(TFwd3D& rData)
   {
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::transferBwd3D(TBwd3D& rData)
   {
      // Convert data
      this->mspConverter3D->convertBwd(rData, this->storage2D());

      // Free input data
      this->storage3D().freeBwd(rData);
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> TFwd3D&  Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::providePhysical()
   {
      return this->storage3D().provideFwd();
   }

   template <typename TFwd1D, typename TBwd1D, typename TFwd2D, typename TBwd2D, typename TFwd3D, typename TBwd3D> void Communicator3D<TFwd1D, TBwd1D, TFwd2D, TBwd2D, TFwd3D, TBwd3D>::holdPhysical(TFwd3D& rData)
   {
      // Hold the input data
      this->storage3D().holdFwd(rData);
   }

}

#endif // COMMUNICATOR3D_HPP
