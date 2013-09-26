/**
 * @file PMIndexConv.hpp
 * @brief Implementation of the index converter for plus-minus FFT frequency ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef PMINDEXCONV_HPP
#define PMINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/TransformResolution.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the index converter for plus-minus FFT frequency ordering
    */
   class PMIndexConv
   {
      public:
         PMIndexConv(SharedCTransformResolution spTRes)
      {mspTRes = spTRes;};
         ~PMIndexConv() {};

         /**
          * @brief Convert first index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert first index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int iS(const int i, const int j, const int k);

         /**
          * @brief Convert first index (2D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int i(const int i, const int j, const int idxI, const int idxJ);

         /**
          * @brief Convert first index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         int iS(const int i, const int j);

         /**
          * @brief Convert first index (1D)
          *
          * @param i    Array index off first dimension
          */
         int i(const int i);

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert second index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int jS(const int i, const int j, const int k);

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         int j(const int i, const int j, const int idxI, const int idxJ);

         /**
          * @brief Convert second index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         int jS(const int i, const int j);

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          * @param idxK Physical index of third dimension
          */
         int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert third index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         int kS(const int i, const int j, const int k);

         SharedCTransformResolution mspTRes;
   };

   inline int PMIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      int nN = this->mspTRes->dim<Dimensions::Data::DAT3D>();

      //std::cerr << this->mspTRes->dim<Dimensions::Data::DAT3D>() << std::endl;

      if (idxK >= nN/2)
      {
         return nN/2 + idxK;
      } else
      {
         return idxK;
      }
   }

   inline int PMIndexConv::iS(const int i, const int j, const int k)
   {
      return PMIndexConv::i(i,j,k,i,j,k);
   }

   inline int PMIndexConv::i(const int i, const int j, const int idxI, const int idxJ)
   {
      return j;
   }

   inline int PMIndexConv::iS(const int i, const int j)
   {
      return PMIndexConv::i(i,j,i,j);
   }

   inline int PMIndexConv::i(const int i)
   {
      return i;
   }

   inline int PMIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return i;
   }

   inline int PMIndexConv::jS(const int i, const int j, const int k)
   {
      return PMIndexConv::j(i,j,k,i,j,k);
   }

   inline int PMIndexConv::j(const int i, const int j, const int idxI, const int idxJ)
   {
      return i;
   }

   inline int PMIndexConv::jS(const int i, const int j)
   {
      return PMIndexConv::j(i,j,i,j);
   }

   inline int PMIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return j;
   }

   inline int PMIndexConv::kS(const int i, const int j, const int k)
   {
      return PMIndexConv::k(i,j,k,i,j,k);
   }

}
}

#endif // PMINDEXCONV_HPP
