/**
 * @file SerialConverterBase.hpp
 * @brief Implementation of the serial data converter base 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SERIALCONVERTERBASE_HPP
#define SERIALCONVERTERBASE_HPP

// Debug includes
//
#include "Profiler/ProfilerMacro.h"
#include "StaticAsserts/StaticAssert.hpp"

// Configuration includes
//

// System includes
//
#include <cassert>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/Dimensions.hpp"
#include "StorageProviders/StoragePairProviderMacro.h"
#include "Communicators/Converters/IConverter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the serial data converter base.
    *
    * \mhdTodo Does index converter really need to be templated ? Is the type selector not enough?
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class SerialConverterBase: public IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>
   {
      public:
         /**
          * @brief Constructor
          */
         SerialConverterBase();

         /**
          * @brief Destructor
          */
         virtual ~SerialConverterBase();

         /**
          * @brief Convert data from TFwdA to TBwdB
          */
         virtual void convertFwd(const TFwdA& in, StoragePairProviderMacro<TFwdB, TBwdB>& storage);

         /**
          * @brief Convert data from TBwdB to TFwdA
          */
         virtual void convertBwd(const TBwdB& in, StoragePairProviderMacro<TFwdA, TBwdA>& storage);
         
      protected:
         /**
          * @brief Get point data from TBwdB scalar (might involve modification of indexes)
          *
          * @param in   Input data
          * @param i    First index of TBwdB extracted from TFwdA
          * @param j    Second index of TBwdB extracted from TFwdA
          * @param k    Third index of TBwdB extracted from TFwdA
          */
         typename TBwdB::PointType bwdPoint(const TBwdB& in, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set point data from TBwdB scalar (might involve modification of indexes)
          *
          * @param rOut Output data
          * @param i    First index of TBwdB extracted from TFwdA
          * @param j    Second index of TBwdB extracted from TFwdA
          * @param k    Third index of TBwdB extracted from TFwdA
          */
         typename TBwdB::PointType& rBwdPoint(TBwdB& rOut, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Store the local transform resolution of the first transform
          */
         SharedCTransformResolution mspTRes;
      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SerialConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::PointType , typename TBwdB::PointType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~SerialConverterBase()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertBwd(const TBwdB& in, StoragePairProviderMacro<TFwdA, TBwdA>& storage)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::BWDSENDCONV);

      // Get storage for output value 
      TFwdA &rOut = storage.provideFwd();

      // 3D case
      if(TFwdB::FieldDimension == Dimensions::THREED)
      {
         // Loop over slowest direction of output
         for(int k = 0; k < this->mspTRes->template dim<Dimensions::Data::DAT3D>(); k++)
         {
            // Loop over slow direction of output
            for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               // Loop over fast direction of output
               for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  rOut.rPoint(i,j,k) = this->bwdPoint(in, i, j, k);
               }
            }
         }

      // 2D case
      } else if(TFwdB::FieldDimension == Dimensions::TWOD)
      {
         // Loop over slow direction of output
         for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(); j++)
         {
            // Loop over fast direction of output
            for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(); i++)
            {
               rOut.rPoint(i,j) = this->bwdPoint(in, i, j);
            }
         }

      // 1D case
      } else if(TFwdB::FieldDimension == Dimensions::ONED)
      {
         //
         // No work is required in 1D
         //
         //
      }

      // Hold the output data
      storage.holdFwd(rOut);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::BWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertFwd(const TFwdA& in, StoragePairProviderMacro<TFwdB, TBwdB>& storage)
   {
      // Start detailed profiler
      DetailedProfilerMacro_start(ProfilerMacro::FWDSENDCONV);

      // Get storage for output value 
      TBwdB &rOut = storage.provideBwd();

      // 3D case
      if(TFwdA::FieldDimension == Dimensions::THREED)
      {
         // Loop over slowest direction of input
         for(int k = 0; k < this->mspTRes->template dim<Dimensions::Data::DAT3D>(); k++)
         {
            // Loop over slow direction of input
            for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(k); j++)
            {
               // Loop over fast direction of input
               for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(k); i++)
               {
                  this->rBwdPoint(rOut, i,j,k) = in.point(i,j,k);
               }
            }
         }

      // 2D case
      } else if(TFwdA::FieldDimension == Dimensions::TWOD)
      {
         // Loop over slow direction of output
         for(int j = 0; j < this->mspTRes->template dim<Dimensions::Data::DAT2D>(); j++)
         {
            for(int i = 0; i < this->mspTRes->template dim<Dimensions::Data::DATF1D>(); i++)
            {
               this->rBwdPoint(rOut, i,j) = in.point(i,j);
            }
         }

      // 1D case
      } else if(TFwdA::FieldDimension == Dimensions::ONED)
      {
         //
         // No work is required in 1D
         //
      }

      // Hold the output data
      storage.holdBwd(rOut);

      // Stop detailed profiler
      DetailedProfilerMacro_stop(ProfilerMacro::FWDSENDCONV);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> typename TBwdB::PointType SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bwdPoint(const TBwdB& in, const int i, const int j, const int k)
   {
      #ifdef QUICC_MPI
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i,k);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j,k);
            int idxK = this->mspTRes->template idx<Dimensions::Data::DAT3D>(k);
            return in.point(this->mspIdxConv->i(i,j,k,idxI,idxJ,idxK),this->mspIdxConv->j(i,j,k,idxI,idxJ,idxK),this->mspIdxConv->k(i,j,k,idxI,idxJ,idxK));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j);
            return in.point(this->mspIdxConv->i(i,j,idxI,idxJ),this->mspIdxConv->j(i,j,idxI,idxJ));

         } else
         {
            return in.point(this->mspIdxConv->i(i));
         }
      #else
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            return in.point(this->mspIdxConv->iS(i,j,k),this->mspIdxConv->jS(i,j,k),this->mspIdxConv->kS(i,j,k));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return in.point(this->mspIdxConv->iS(i,j),this->mspIdxConv->jS(i,j));

         } else
         {
            return in.point(this->mspIdxConv->i(i));
         }
      #endif //QUICC_MPI
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> typename TBwdB::PointType& SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::rBwdPoint(TBwdB& rOut, const int i, const int j, const int k)
   {
      #ifdef QUICC_MPI
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i,k);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j,k);
            int idxK = this->mspTRes->template idx<Dimensions::Data::DAT3D>(k);

            return rOut.rPoint(this->mspIdxConv->i(i,j,k,idxI,idxJ,idxK),this->mspIdxConv->j(i,j,k,idxI,idxJ,idxK),this->mspIdxConv->k(i,j,k,idxI,idxJ,idxK));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            int idxI = this->mspTRes->template idx<Dimensions::Data::DATF1D>(i);
            int idxJ = this->mspTRes->template idx<Dimensions::Data::DAT2D>(j);
            return rOut.rPoint(this->mspIdxConv->i(i,j,idxI,idxJ),this->mspIdxConv->j(i,j,idxI,idxJ));

         } else
         {
            return rOut.rPoint(this->mspIdxConv->i(i));
         }
      #else
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            return rOut.rPoint(this->mspIdxConv->iS(i,j,k),this->mspIdxConv->jS(i,j,k),this->mspIdxConv->kS(i,j,k));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return rOut.rPoint(this->mspIdxConv->iS(i,j),this->mspIdxConv->jS(i,j));

         } else
         {
            return rOut.rPoint(this->mspIdxConv->i(i));
         }
      #endif //QUICC_MPI
   }

}
}

#endif // SERIALCONVERTERBASE_HPP
