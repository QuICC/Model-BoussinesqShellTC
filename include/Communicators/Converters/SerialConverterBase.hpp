/** \file SerialConverterBase.hpp
 *  \brief Implementation of the serial data converter base
 *
 *  \mhdBug Needs test
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
#include "StorageProviders/StoragePairProvider.hpp"
#include "Communicators/Converters/IConverter.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * \brief Implementation of the serial data converter base.
    */
   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> class SerialConverterBase: public IConverter<TFwdA, TBwdA, TFwdB, TBwdB>
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
         virtual void convertFwd(const TFwdA& in, StoragePairProvider<TFwdB, TBwdB>& storage);

         /**
          * @brief Convert data from TBwdB to TFwdA
          */
         virtual void convertBwd(const TBwdB& in, StoragePairProvider<TFwdA, TBwdA>& storage);
         
      protected:
         /**
          * @brief Get point data from TBwdB scalar (might involve modification of indexes)
          *
          * @param in   Input data
          * @param i    First index of TBwdB extracted from TFwdA
          * @param j    Second index of TBwdB extracted from TFwdA
          * @param k    Third index of TBwdB extracted from TFwdA
          */
         const typename TBwdB::CoefficientType& bwdPoint(const TBwdB& in, const int i, const int j = 0, const int k = 0);

         /**
          * @brief Set point data from TBwdB scalar (might involve modification of indexes)
          *
          * @param rOut Output data
          * @param i    First index of TBwdB extracted from TFwdA
          * @param j    Second index of TBwdB extracted from TFwdA
          * @param k    Third index of TBwdB extracted from TFwdA
          */
         typename TBwdB::CoefficientType& rBwdPoint(TBwdB& rOut, const int i, const int j = 0, const int k = 0);

         SharedTransformResolution mspTRes;

      private:
   };

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::SerialConverterBase()
      : IConverter<TFwdA, TBwdA, TFwdB, TBwdB>()
   {
      // Check that all dimensions match
      Debug::StaticAssert< (TFwdA::FieldDimension == TBwdA::FieldDimension) >();
      Debug::StaticAssert< (TBwdA::FieldDimension == TFwdB::FieldDimension) >();
      Debug::StaticAssert< (TFwdB::FieldDimension == TBwdB::FieldDimension) >();

      // Check that the data type is the same
      Debug::StaticTypeAssert<typename TFwdA::CoefficientType , typename TBwdB::CoefficientType>();
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::~SerialConverterBase()
   {
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertBwd(const TBwdB& in, StoragePairProvider<TFwdA, TBwdA>& storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::BWDCONVSEND);

      // Get storage for output value 
      TFwdA &rOut = storage.provideFwd();

      // 3D case
      if(TFwdB::FieldDimension == Dimensions::THREED)
      {
         // Loop over slowest direction of output
         for(int k = 0; k < rOut.dim3D(); k++)
         {
            // Loop over slow direction of output
            for(int j = 0; j < rOut.dim2D(k); j++)
            {
               // Loop over fast direction of output
               for(int i = 0; i < rOut.dim1D(j,k); i++)
               {
                  rOut.rPoint(i,j,k) = this->bwdPoint(in, i, j, k);
               }
            }
         }

      // 2D case
      } else if(TFwdB::FieldDimension == Dimensions::TWOD)
      {
         // Loop over slow direction of output
         for(int j = 0; j < rOut.dim2D(); j++)
         {
            // Loop over fast direction of output
            for(int i = 0; i < rOut.dim1D(j); i++)
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

      DetailedProfilerMacro_stop(ProfilerMacro::BWDCONVSEND);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> void SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::convertFwd(const TFwdA& in, StoragePairProvider<TFwdB, TBwdB>& storage)
   {
      DetailedProfilerMacro_start(ProfilerMacro::FWDCONVSEND);

      // Get storage for output value 
      TBwdB &rOut = storage.provideBwd();

      // 3D case
      if(TFwdA::FieldDimension == Dimensions::THREED)
      {
         // Loop over slowest direction of input
         for(int k = 0; k < in.dim3D(); k++)
         {
            // Loop over slow direction of input
            for(int j = 0; j < in.dim2D(k); j++)
            {
               // Loop over fast direction of input
               for(int i = 0; i < in.dim1D(j,k); i++)
               {
                  this->rBwdPoint(rOut, i,j,k) = in.point(i,j,k);
               }
            }
         }

      // 2D case
      } else if(TFwdA::FieldDimension == Dimensions::TWOD)
      {
         // Loop over slow direction of output
         for(int j = 0; j < in.dim2D(); j++)
         {
            for(int i = 0; i < in.dim1D(j); i++)
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

      DetailedProfilerMacro_stop(ProfilerMacro::FWDCONVSEND);
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> const typename TBwdB::CoefficientType& SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::bwdPoint(const TBwdB& in, const int i, const int j, const int k)
   {
      #ifdef GEOMHDISCC_MPI
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            int idxI = this->mspTRes->idxFwd(i,j,k);
            int idxJ = this->mspTRes->idx2D(j,k);
            int idxK = this->mspTRes->idx3D(k);
            return in.point(TIdx::i(i,j,k,idxI,idxJ,idxK),TIdx::j(i,j,k,idxI,idxJ,idxK),TIdx::k(i,j,k,idxI,idxJ,idxK));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            int idxI = this->mspTRes->idxFwd(i,j);
            int idxJ = this->mspTRes->idx2D(j);
            return in.point(TIdx::i(i,j,idxI,idxJ),TIdx::j(i,j,idxI,idxJ));

         } else if(TBwdB::FieldDimension == Dimensions::ONED)
         {
            return in.point(TIdx::i(i));
         }
      #else
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            return in.point(TIdx::iS(i,j,k),TIdx::jS(i,j,k),TIdx::kS(i,j,k));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return in.point(TIdx::iS(i,j),TIdx::jS(i,j));

         } else if(TBwdB::FieldDimension == Dimensions::ONED)
         {
            return in.point(TIdx::i(i));
         }
      #endif //GEOMHDISCC_MPI
   }

   template <typename TFwdA, typename TBwdA, typename TFwdB, typename TBwdB, typename TIdx> typename TBwdB::CoefficientType& SerialConverterBase<TFwdA, TBwdA, TFwdB, TBwdB, TIdx>::rBwdPoint(TBwdB& rOut, const int i, const int j, const int k)
   {
      #ifdef GEOMHDISCC_MPI
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            int idxI = this->mspTRes->idxFwd(i,j,k);
            int idxJ = this->mspTRes->idx2D(j,k);
            int idxK = this->mspTRes->idx3D(k);

            return rOut.rPoint(TIdx::i(i,j,k,idxI,idxJ,idxK),TIdx::j(i,j,k,idxI,idxJ,idxK),TIdx::k(i,j,k,idxI,idxJ,idxK));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            int idxI = this->mspTRes->idxFwd(i,j);
            int idxJ = this->mspTRes->idx2D(j);
            return rOut.rPoint(TIdx::i(i,j,idxI,idxJ),TIdx::j(i,j,idxI,idxJ));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return rOut.rPoint(TIdx::i(i));
         }
      #else
         if(TBwdB::FieldDimension == Dimensions::THREED)
         {
            return rOut.rPoint(TIdx::iS(i,j,k),TIdx::jS(i,j,k),TIdx::kS(i,j,k));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return rOut.rPoint(TIdx::iS(i,j),TIdx::jS(i,j));

         } else if(TBwdB::FieldDimension == Dimensions::TWOD)
         {
            return rOut.rPoint(TIdx::i(i));
         }
      #endif //GEOMHDISCC_MPI
   }

}
}

#endif // SERIALCONVERTERBASE_HPP
