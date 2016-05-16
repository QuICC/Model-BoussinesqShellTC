/** 
 * @file WorlandPolynomial.hpp
 * @brief Implementation of the Worland polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef WORLANDPOLYNOMIAL_HPP
#define WORLANDPOLYNOMIAL_HPP

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
    * @brief Implementation of the Worland polynomial
    */ 
   class WorlandPolynomial
   {
      public:
         /**
          * @brief Compute \f$W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Wnl(Matrix& poly, internal::Matrix& ipoly, const int l, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d}{dr} W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{W_n^l (r)}{r}\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void Wnl_r(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

      private:
         /**
          * @brief Get alpha parameter of Jacobi polynomial
          */
         static internal::MHDFloat alpha(const int l);

         /**
          * @brief Get beta parameter of Jacobi polynomial
          */
         static internal::MHDFloat beta(const int l);

         /**
          * @brief Multiply Jacobi polynomial by \f$r^l\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void rl(internal::Matrix& ipoly, const int l, const internal::Array& igrid);

         /**
          * @brief Constructor
          */
         WorlandPolynomial();

         /**
          * @brief Destructor
          */
         ~WorlandPolynomial();

   };
}
}

#endif // WORLANDPOLYNOMIAL_HPP
