/**
 * @file NoIndexConv.hpp
 * @brief Implementation of the index converter doing nothing 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef NOINDEXCONV_HPP
#define NOINDEXCONV_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//

namespace GeoMHDiSCC {

namespace Parallel {

   /**
    * \brief Implementation of the index converter doing nothing
    */
   class NoIndexConv
   {
      public:
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
         static int i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert first index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         static int iS(const int i, const int j, const int k);

         /**
          * @brief Convert first index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         static int i(const int i, const int j, const int idxI, const int idxJ);

         /**
          * @brief Convert first index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         static int iS(const int i, const int j);

         /**
          * @brief Convert first index (1D)
          *
          * @param i    Array index off first dimension
          */
         static int i(const int i);

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
         static int j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert second index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         static int jS(const int i, const int j, const int k);

         /**
          * @brief Convert second index (3D)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param idxI Physical index of first dimension
          * @param idxJ Physical index of second dimension
          */
         static int j(const int i, const int j, const int idxI, const int idxJ);

         /**
          * @brief Convert second index (2D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          */
         static int jS(const int i, const int j);

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
         static int k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK);

         /**
          * @brief Convert third index (3D) (specialised serial version)
          *
          * @param i    Array index off first dimension
          * @param j    Array index of second dimension
          * @param k    Array index of third dimension
          */
         static int kS(const int i, const int j, const int k);
   };

   inline int NoIndexConv::i(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return k;
   }

   inline int NoIndexConv::iS(const int i, const int j, const int k)
   {
      return NoIndexConv::i(i,j,k,i,j,k);
   }

   inline int NoIndexConv::i(const int i, const int j, const int idxI, const int idxJ)
   {
      return j;
   }

   inline int NoIndexConv::iS(const int i, const int j)
   {
      return NoIndexConv::i(i,j,i,j);
   }

   inline int NoIndexConv::i(const int i)
   {
      return i;
   }

   inline int NoIndexConv::j(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return i;
   }

   inline int NoIndexConv::jS(const int i, const int j, const int k)
   {
      return NoIndexConv::j(i,j,k,i,j,k);
   }

   inline int NoIndexConv::j(const int i, const int j, const int idxI, const int idxJ)
   {
      return i;
   }

   inline int NoIndexConv::jS(const int i, const int j)
   {
      return NoIndexConv::j(i,j,i,j);
   }

   inline int NoIndexConv::k(const int i, const int j, const int k, const int idxI, const int idxJ, const int idxK)
   {
      return j;
   }

   inline int NoIndexConv::kS(const int i, const int j, const int k)
   {
      return NoIndexConv::k(i,j,k,i,j,k);
   }

}
}

#endif // NOINDEXCONV_HPP
