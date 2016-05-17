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

         enum NormalizationId {
            NATURAL,
            UNITWORLAND,
         };

         /**
          * @brief Set the normalization
          */
         static void setNormalization(const NormalizationId id);

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
          * @brief Compute \f$e P_n^{(\alpha,\beta)} (x) + 2(x+1)\frac{d}{dr}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drrPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat e, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d}{dr}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

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
          * @brief Compute \f$c P_n^{(\alpha,\beta)}(x) + 2(x+1)\frac{d}{dx}P_n^{(\alpha,\beta)}(x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drrPnab(Eigen::Ref<internal::Matrix> idrrpnab, const internal::MHDFloat e, const Eigen::Ref<const internal::Matrix>& ipnab, const Eigen::Ref<const internal::Matrix>& idpnab, const internal::Array& igrid);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drrP0ab(Eigen::Ref<internal::Matrix> idrrpnab, const internal::MHDFloat e, const Eigen::Ref<const internal::Matrix>& ipnab);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drPnab(Eigen::Ref<internal::Matrix> idrpnab, const Eigen::Ref<const internal::Matrix>& idnab, const internal::Array& igrid);

         /**
          * @brief Polynomial normalizer for natural normalization
          */
         static internal::Array naturalPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 0th polynomial normalizer for natural normalization
          */
         static internal::Array naturalP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 1st polynomial normalizer for natural normalization
          */
         static internal::Array naturalP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative normalizer for natural normalization
          */
         static internal::Array naturalDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 0th derivative normalizer for natural normalization
          */
         static internal::Array naturalDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 1st first derivative normalizer for natural normalization
          */
         static internal::Array naturalDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Polynomial normalizer for unit Worland normalization
          */
         static internal::Array unitWPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 0th polynomial normalizer for unit Worland normalization
          */
         static internal::Array unitWP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 1st polynomial normalizer for unit Worland normalization
          */
         static internal::Array unitWP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 0th first derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief 1st first derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         static NormalizerNAB mpNormPnab;
         static NormalizerAB mpNormP0ab;
         static NormalizerAB mpNormP1ab;

         static NormalizerNAB mpNormDPnab;
         static NormalizerAB mpNormDP0ab;
         static NormalizerAB mpNormDP1ab;

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
