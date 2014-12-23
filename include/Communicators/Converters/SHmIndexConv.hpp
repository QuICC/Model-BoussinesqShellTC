/**
 * @file SHmIndexConv.hpp
 * @brief Implementation of the index converter for spherical harmonics with m spectral ordering
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHMINDEXCONV_HPP
#define SHMINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * @brief Implementation of the index converter for spherical harmonics with m spectral ordering
    */
   class SHmIndexConv
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes   Shared resolution
          * @param id      Forward dimension index ID
          */
         SHmIndexConv(SharedResolution spRes, const Dimensions::Transform::Id id);

         /**
          * @brief Destructor
          */
         ~SHmIndexConv();

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

   inline int SHmIndexConv::centralPadding(const int idx, const int k)
   {
      return 0;
   }

   inline int SHmIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return idxK+idxJ;
   }

   inline int SHmIndexConv::iS(const int i, const int j, const int k)
   {
      return SHmIndexConv::i(i,j,k,i,j,k);
   }

   inline int SHmIndexConv::i(const int i, const int j, const int idxI, const int idxJ)
   {
      return idxJ+idxI;
   }

   inline int SHmIndexConv::iS(const int i, const int j)
   {
      return SHmIndexConv::i(i,j,i,j);
   }

   inline int SHmIndexConv::i(const int i)
   {
      return i;
   }

   inline int SHmIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return i;
   }

   inline int SHmIndexConv::jS(const int i, const int j, const int k)
   {
      return SHmIndexConv::j(i,j,k,i,j,k);
   }

   inline int SHmIndexConv::j(const int i, const int j, const int idxI, const int idxJ)
   {
      return i;
   }

   inline int SHmIndexConv::jS(const int i, const int j)
   {
      return SHmIndexConv::j(i,j,i,j);
   }

   inline int SHmIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return k;
   }

   inline int SHmIndexConv::kS(const int i, const int j, const int k)
   {
      return SHmIndexConv::k(i,j,k,i,j,k);
   }

}
}

#endif // SHMINDEXCONV_HPP
