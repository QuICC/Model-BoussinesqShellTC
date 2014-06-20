/**
 * @file SHIndexConv.hpp
 * @brief Implementation of the index converter for spherical harmonics 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHINDEXCONV_HPP
#define SHINDEXCONV_HPP

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
    * @brief Implementation of the index converter for spherical harmonics
    */
   class SHIndexConv
   {
      public:
         /**
          * @brief Constructor
          *
          * @param spRes   Shared resolution
          * @param id      Forward dimension index ID
          */
         SHIndexConv(SharedResolution spRes, const Dimensions::Transform::Id id);

         /**
          * @brief Destructor
          */
         ~SHIndexConv();

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

   inline int SHIndexConv::centralPadding(const int idx, const int k)
   {
      return 0;
   }

   inline int SHIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return idxK-idxJ;
   }

   inline int SHIndexConv::iS(const int i, const int j, const int k)
   {
      return SHIndexConv::i(i,j,k,i,j,k);
   }

   inline int SHIndexConv::i(const int i, const int j, const int idxI, const int idxJ)
   {
      return idxJ-idxI;
   }

   inline int SHIndexConv::iS(const int i, const int j)
   {
      return SHIndexConv::i(i,j,i,j);
   }

   inline int SHIndexConv::i(const int i)
   {
      return i;
   }

   inline int SHIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return i;
   }

   inline int SHIndexConv::jS(const int i, const int j, const int k)
   {
      return SHIndexConv::j(i,j,k,i,j,k);
   }

   inline int SHIndexConv::j(const int i, const int j, const int idxI, const int idxJ)
   {
      return i;
   }

   inline int SHIndexConv::jS(const int i, const int j)
   {
      return SHIndexConv::j(i,j,i,j);
   }

   inline int SHIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return j;
   }

   inline int SHIndexConv::kS(const int i, const int j, const int k)
   {
      return SHIndexConv::k(i,j,k,i,j,k);
   }

}
}

#endif // SHINDEXCONV_HPP
