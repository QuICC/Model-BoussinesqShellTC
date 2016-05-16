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
         static void Pnab(Eigen::Ref<internal::Matrix> iplm, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid);

         /**
          * @brief Compute \f$P_0^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P0ab(Eigen::Ref<internal::Matrix> ip0ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$P_1^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void P1ab(Eigen::Ref<internal::Matrix> ip1ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ip0ab, const internal::Array& igrid);

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
