/**
 * @file DecoupledComplexInternal.hpp
 * @brief Useful methods for the DecoupledComplex type
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DECOUPLEDCOMPLEXINTERNAL_HPP
#define DECOUPLEDCOMPLEXINTERNAL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"

namespace GeoMHDiSCC {

namespace Datatypes  {

   namespace internal
   {
      /**
       * @brief Get stored value as complex number from DecoupledComplex storage
       *
       * @param mat
       * @param k
       */
      MHDComplex getScalar(const DecoupledZMatrix& mat, const int k);

      /**
       * @brief Get stored value 
       *
       * @tparam TData
       * @param mat
       * @param k
       */
      template <typename TData> typename TData::Scalar getScalar(const TData& mat, const int k);

      /**
       * @brief Set complex number into DecoupledComplex storage 
       *
       * @param mat
       * @param k
       * @param val
       */
      void setScalar(DecoupledZMatrix& mat, const int k, const MHDComplex& val);

      /**
       * @brief Set value 
       *
       * @tparam TData
       * @param mat
       * @param k
       * @param val
       */
      template <typename TData> void setScalar(TData& mat, const int k, const typename TData::Scalar& val);

      /**
       * @brief Add complex value to DecoupledComplex storage 
       *
       * @param mat
       * @param k
       * @param val
       */
      void addScalar(DecoupledZMatrix& mat, const int k, const MHDComplex& val);

      /**
       * @brief Add value
       *
       * @tparam TData
       * @param mat
       * @param k
       * @param val
       */
      template <typename TData> void addScalar(TData& mat, const int k, const typename TData::Scalar& val);


      inline MHDComplex getScalar(const DecoupledZMatrix& mat, const int k)
      {
         return MHDComplex(mat.real()(k),mat.imag()(k));
      }

      template <typename TData> inline typename TData::Scalar getScalar(const TData& mat, const int k)
      {
         return mat(k);
      }

      inline void setScalar(DecoupledZMatrix& mat, const int k, const MHDComplex& val)
      {
         mat.real()(k) = val.real();
         mat.imag()(k) = val.imag();
      }

      template <typename TData> inline void setScalar(TData& mat, const int k, const typename TData::Scalar& val)
      {
         mat(k) = val;
      }

      inline void addScalar(DecoupledZMatrix& mat, const int k, const MHDComplex& val)
      {
         mat.real()(k) += val.real();
         mat.imag()(k) += val.imag();
      }

      template <typename TData> inline void addScalar(TData& mat, const int k, const typename TData::Scalar& val)
      {
         mat(k) += val;
      }
   }
}
}

#endif // DECOUPLEDCOMPLEXINTERNAL_HPP
