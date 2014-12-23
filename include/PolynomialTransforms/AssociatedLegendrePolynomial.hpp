/** 
 * @file AssociatedLegendrePolynomial.hpp
 * @brief Implementation of the associated Legendre polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef ASSOCIATEDLEGENDREPOLYNOMIAL_HPP
#define ASSOCIATEDLEGENDREPOLYNOMIAL_HPP

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
    * @brief Implementation of the associated Legendre polynomial
    */ 
   class AssociatedLegendrePolynomial
   {
      public:
         /**
          * @brief Constructor
          */
         AssociatedLegendrePolynomial();

         /**
          * @brief Destructor
          */
         ~AssociatedLegendrePolynomial();

         /**
          * @brief Compute the associated Legendre \f$P_m^l (\cos\theta)\f$ for all l
          *
          * Internal computation can be done in multiple precision
          */
         static void Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_m^l (\cos\theta)\f$ for all l
          *
          * Internal computation can be done in multiple precision
          */
         static void dPlm(Matrix& diff, internal::Matrix& idiff, const int m, const internal::Matrix& ipoly, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{P_m^l (\cos\theta)}{\sin\theta}\f$ for all l
          *
          * Internal computation can be done in multiple precision
          */
         static void sin_1Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$P_m^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pmm(internal::Array& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$P_m^{m+1} (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_m^m (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm(internal::Array& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_\theta} P_m^{m+1} (\cos\theta)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& idpmm, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{P_m^m(\cos\theta)}{\sin\theta}\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void sin_1Pmm(internal::Array& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{P_m^{m+1}(\cos\theta)}{\sin\theta}\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void sin_1Pmm1(internal::Array& op, const int m, const internal::Array& isin_1pmm, const internal::Array& igrid);

      private:
   };
}
}

#endif // ASSOCIATEDLEGENDREPOLYNOMIAL_HPP
