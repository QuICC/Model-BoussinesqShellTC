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
          * @brief Compute the associated Legendre \f$P_m^l (x)\f$ for all l
          *
          * Internal computation can be done in multiple precision
          */
         static void Plm(Matrix& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_x} P_m^l (x)\f$ for all l
          *
          * Internal computation can be done in multiple precision
          */
         static void dPlm(Matrix& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$P_m^m (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pmm(internal::Array& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$P_m^{m+1} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Pmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_x} P_m^m (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm(internal::Array& op, const int m, const internal::Array& igrid);

         /**
          * @brief Compute the associated Legendre \f$\frac{d}{d_x} P_m^{m+1} (x)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dPmm1(internal::Array& op, const int m, const internal::Array& igrid);

      private:
   };
}
}

#endif // ASSOCIATEDLEGENDREPOLYNOMIAL_HPP
