/**
 * @file MpiConverterBase.hpp
 * @brief Templated implementation of the base of a MPI data converter 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef MPICONVERTERBASE_HPP
#define MPICONVERTERBASE_HPP

// Debug includes
//
#include "StorageProfiler/StorageProfilerMacro.h"
#include "Profiler/ProfilerMacro.h"

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <cassert>
#include <set>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MpiTypes.hpp"
#include "Enums/TransformDirection.hpp"
#include "Resolutions/Resolution.hpp"
#include "Communicators/Converters/IConverter.hpp"
#include "Communicators/Converters/MpiConverterTools.hpp"
#include "Communicators/CommunicationBuffer.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Templated implementation of the base of a MPI data converter.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class MpiConverterBase: public IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         #if defined QUICC_MPIPACK_MANUAL
            typedef SharedPtrMacro<CommunicationBuffer<typename TFwdA::PointType> > SharedFwdBufferType;
            typedef SharedPtrMacro<CommunicationBuffer<typename TBwdB::PointType> > SharedBwdBufferType;
         #else
            typedef SharedPtrMacro<CommunicationBuffer<char> > SharedFwdBufferType;
            typedef SharedPtrMacro<CommunicationBuffer<char> > SharedBwdBufferType;
         #endif //defined QUICC_MPIPACK_MANUAL

         /**
          * @brief Constructor
          */
         MpiConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~MpiConverterBase();

         /**
          * @brief Set the communication buffers
          *
          * @brief spFwd Forward communication buffers
          * @brief spBwd Backward communication buffers
          */
         void setBuffers(SharedFwdBufferType spFwd, SharedBwdBufferType spBwd);

         /**
          * @brief Get forward buffer sizes
          */
         const std::vector<int> & fwdSizes() const;

         /**
          * @brief Get backward buffer sizes
          */
         const std::vector<int> & bwdSizes() const;
         
      protected:
         /**
          * @brief Reset Fwd buffer positions
          */
         void resetFwdPositions();

         /**
          * @brief Reset Bwd buffer positions
          */
         void resetBwdPositions();

         /**
          * @brief Size of the forward packet
          *
          * @param id ID of the node
          */
         int sizeFPacket(const int id) const;

         /**
          * @brief Size of the backward packet
          *
          * @param id ID of the node
          */
         int sizeBPacket(const int id) const;

         /**
          * @brief Get size of the forward CPU group
          */
         int nFCpu() const;

         /**
          * @brief Get size of the backward CPU group
          */
         int nBCpu() const;

         /**
          * @brief Get MPI rank of CPU from forward CPU group
          *
          * @param id CPU group id
          */
         int fCpu(const int id) const;

         /**
          * @brief Get MPI rank of CPU from backward CPU group
          *
          * @param id CPU group id
          */
         int bCpu(const int id) const;

         /**
          * @brief Keep empty communcations?
          */
         bool mNeedEmptyComm;

         /**
          * @brief Communication packs counter
          */
         int mPacks;

         /**
          * @brief Direction of operations
          */
         TransformDirection::Id   mDirection;

         /**
          * @brief Transform ID
          */
         Dimensions::Transform::Id mTraId;

         /**
          * @brief List of CPU ranks involved in the forward conversion
          */
         std::vector<int>  mFCpuGroup;

         /**
          * @brief List of CPU ranks involved in the backward conversion
          */
         std::vector<int>  mBCpuGroup;

         /**
          * @brief Forward communication
          */
         SharedFwdBufferType mspFBuffers;

         /**
          * @brief Backward communication buffer pointer
          */
         SharedBwdBufferType mspBBuffers;

         /**
          * @brief List of the forward buffer sizes
          */
         std::vector<int>  mFSizes;

         /**
          * @brief List of the backward buffer sizes
          */
         std::vector<int>  mBSizes;

         /**
          * @brief Possible forward transform packs
          */
         ArrayI   mForwardPacks;

         /**
          * @brief Possible backward transform packs
          */
         ArrayI   mBackwardPacks;

         #if defined QUICC_MPIPACK_MANUAL
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<std::vector<typename MpiConverterTools<TFwdA::FieldDimension>::Coordinate> >  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<std::vector<typename MpiConverterTools<TBwdB::FieldDimension>::Coordinate> > mBTypes;

         #else
            /**
             * @brief Storage for the forward datatypes
             */
            std::vector<MPI_Datatype>  mFTypes;

            /**
             * @brief Storage for the backward datatypes
             */
            std::vector<MPI_Datatype> mBTypes;
         #endif //defined QUICC_MPIPACK_MANUAL

      private:
         /**
          * @brief Cleanup the data types
          */
         void cleanupTypes();
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline const std::vector<int>& MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::fwdSizes() const
   {
      return this->mFSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline const std::vector<int>& MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bwdSizes() const
   {
      return this->mBSizes;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::setBuffers(typename MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SharedFwdBufferType spFwd, typename MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SharedBwdBufferType spBwd)
   {
      // Set the forward buffers
      this->mspFBuffers = spFwd;

      // Set the backward buffers
      this->mspBBuffers = spBwd;
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::nFCpu() const
   {
      return this->mFCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::nBCpu() const
   {
      return this->mBCpuGroup.size();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::fCpu(const int id) const
   {
      return this->mFCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bCpu(const int id) const
   {
      return this->mBCpuGroup.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sizeFPacket(const int id) const
   {
      return this->mPacks*this->mFSizes.at(id);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> inline int MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::sizeBPacket(const int id) const
   {
      return this->mPacks*this->mBSizes.at(id);
   }
      
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::MpiConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>(), mNeedEmptyComm(false), mPacks(0)
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~MpiConverterBase()
   {
      // Cleanup mpi datatypes 
      void cleanupTypes();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::resetFwdPositions()
   {
      this->mspFBuffers->resetPositions();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::resetBwdPositions()
   {
      this->mspBBuffers->resetPositions();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void MpiConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::cleanupTypes()
   {
      #if defined QUICC_MPIPACK_MANUAL
         // No cleanup is required
         
      #else
         // Cleanup the F Types
         typename std::vector<MPI_Datatype>::iterator  it;

         for(it = this->mFTypes.begin(); it != this->mFTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }

         for(it = this->mBTypes.begin(); it != this->mBTypes.end(); ++it)
         {
            MPI_Type_free(&(*it));
         }
      #endif //defined QUICC_MPIPACK_MANUAL
   }

}
}

#endif // MPICONVERTERBASE_HPP
