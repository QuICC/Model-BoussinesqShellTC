/** 
 * @file JacobiPolynomial.hpp
 * @brief Implementation of the Jacobi polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef JACOBIPOLYNOMIAL_HPP
#define JACOBIPOLYNOMIAL_HPP

// Debug includes
//

// Configuration includes
//

// System includes
//
#include <tr1/functional>

// External includes
//

// Project includes
//
#include "Base/Precision.hpp"

namespace GeoMHDiSCC {

namespace Polynomial {

   /**
    * @brief Implementation of the Jacobi polynomial
    */ 
   class JacobiPolynomial
   {
      public:
         /// Typedef for the function signature of an n independent normalizer 
         typedef std::tr1::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat)> NormalizerAB;
         /// Typedef for the function signature of an n dependent normalizer 
         typedef std::tr1::function<internal::Array(const internal::MHDFloat, const internal::MHDFloat, const internal::MHDFloat)> NormalizerNAB;

         /**
          * @brief Compute \f$P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pnab(Matrix& poly, internal::Matrix& ipoly, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d}{dx}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d^2}{dx^2}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void d2Pnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pnab(Eigen::Ref<internal::Matrix> iplm, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid, NormalizerNAB norm);

         /**
          * @brief Compute \f$P_0^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P0ab(Eigen::Ref<internal::Matrix> ip0ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1ab(Eigen::Ref<internal::Matrix> ip1ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ip0ab, const internal::Array& igrid, NormalizerAB norm);

         /**
          * @brief Polynomial normalizer for natural normalization
          */
         static internal::Array naturalPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Polynomial n=0 normalizer for natural normalization
          */
         static internal::Array naturalP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Polynomial n=1 normalizer for natural normalization
          */
         static internal::Array naturalP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative normalizer for natural normalization
          */
         static internal::Array naturalDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative normalizer for natural normalization
          */
         static internal::Array naturalD2Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalD2P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalD2P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative normalizer for natural normalization
          */
         static internal::Array naturalD3Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalD3P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalD3P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

      private:
         /**
          * @brief Constructor
          */
         JacobiPolynomial();

         /**
          * @brief Destructor
          */
         ~JacobiPolynomial();

   };
}
}

#endif // JACOBIPOLYNOMIAL_HPP
