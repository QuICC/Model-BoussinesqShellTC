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
#include "PolynomialTransforms/ThreeTermRecurrence.hpp"

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
          * @brief Compute \f$\frac{d^2}{dx^2}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void d2Pnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d^3}{dx^3}P_n^{(\alpha,\beta)} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void d3Pnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid);
         static ThreeTermRecurrence::NormalizerNAB normPnab();
         static ThreeTermRecurrence::NormalizerAB normP1ab();
         static ThreeTermRecurrence::NormalizerAB normP0ab();
         static ThreeTermRecurrence::NormalizerNAB normDPnab();
         static ThreeTermRecurrence::NormalizerAB normDP1ab();
         static ThreeTermRecurrence::NormalizerAB normDP0ab();
         static ThreeTermRecurrence::NormalizerNAB normD2Pnab();
         static ThreeTermRecurrence::NormalizerAB normD2P1ab();
         static ThreeTermRecurrence::NormalizerAB normD2P0ab();
         static ThreeTermRecurrence::NormalizerNAB normD3Pnab();
         static ThreeTermRecurrence::NormalizerAB normD3P1ab();
         static ThreeTermRecurrence::NormalizerAB normD3P0ab();

         /**
          * @brief Polynomial normalizer for natural normalization
          */
         static internal::Array naturalPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=0 normalizer for natural normalization
          */
         static internal::Array naturalP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Polynomial n=1 normalizer for natural normalization
          */
         static internal::Array naturalP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative normalizer for natural normalization
          */
         static internal::Array naturalDPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalDP0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief First derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalDP1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative normalizer for natural normalization
          */
         static internal::Array naturalD2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Second derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative normalizer for natural normalization
          */
         static internal::Array naturalD3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=0 normalizer for natural normalization
          */
         static internal::Array naturalD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b);

         /**
          * @brief Third derivative n=1 normalizer for natural normalization
          */
         static internal::Array naturalD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b);

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
