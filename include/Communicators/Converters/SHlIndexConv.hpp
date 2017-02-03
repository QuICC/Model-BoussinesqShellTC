/**
 * @file SHlIndexConv.hpp
 * @brief Implementation of the index converter for spherical harmonics with l spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHLINDEXCONV_HPP
#define SHLINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"

namespace QuICC {

namespace Parallel {

   /**
    * @brief Implementation of the index converter for spherical harmonics with l spectral ordering
    */
   class SHlIndexConv
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes   Shared resolution
          * @param id      Forward dimension index ID
          */
         SHlIndexConv(SharedResolution spRes, const Dimensions::Transform::Id id);

         /**
          * @brief Destructor
          */
         ~SHlIndexConv();

         /**
          * @brief Compute shift due to central padding
          *
          * @param idx  Reference index
          * @param k    Array index of third dimension
          */
         int centralPadding(const int idx, const int k);

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
          * @brief Convert second index (2D)
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
   };

   inline int SHlIndexConv::centralPadding(const int idx, const int k)
   {
      return 0;
   }

   inline int SHlIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return idxK - idxJ;
   }

   inline int SHlIndexConv::iS(const int i, const int j, const int k)
   {
      return SHlIndexConv::i(i,j,k,i,j,k);
   }

   inline int SHlIndexConv::i(const int i, const int j, const int idxI, const int idxJ)
   {
      return idxJ - idxI;
   }

   inline int SHlIndexConv::iS(const int i, const int j)
   {
      return SHlIndexConv::i(i,j,i,j);
   }

   inline int SHlIndexConv::i(const int i)
   {
      return i;
   }

   inline int SHlIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return i;
   }

   inline int SHlIndexConv::jS(const int i, const int j, const int k)
   {
      return SHlIndexConv::j(i,j,k,i,j,k);
   }

   inline int SHlIndexConv::j(const int i, const int j, const int idxI, const int idxJ)
   {
      return i;
   }

   inline int SHlIndexConv::jS(const int i, const int j)
   {
      return SHlIndexConv::j(i,j,i,j);
   }

   inline int SHlIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return j;
   }

   inline int SHlIndexConv::kS(const int i, const int j, const int k)
   {
      return SHlIndexConv::k(i,j,k,i,j,k);
   }

}
}

#endif // SHLINDEXCONV_HPP
