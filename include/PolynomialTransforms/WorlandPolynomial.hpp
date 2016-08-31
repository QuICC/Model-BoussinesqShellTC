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
#include "PolynomialTransforms/JacobiPolynomial.hpp"

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
          * @brief Compute \f$\frac{d}{dr} r W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void drWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{1}{r}\frac{d}{dr} r W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void r_1drWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{d}{dr} W_n^0 (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dWn0(Matrix& diff, internal::Matrix& idiff, const internal::Array& igrid);

         /**
          * @brief Compute \f$\frac{W_n^l (r)}{r}\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void r_1Wnl(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

         /**
          * @brief Compute spherical \f$\nabla^2 W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void slaplWnl(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

         /**
          * @brief Compute cylindrical horizontal \f$\nabla_h^2 W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void claplhWnl(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

         /**
          * @brief Compute cylindrical horizontal \f$\frac{1}{r}\nabla_h^2 W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void r_1claplhWnl(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

         /**
          * @brief Compute cylindrical horizontal \f$\frac{d}{dr}\nabla_h^2 W_n^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void dclaplhWnl(Matrix& poly, internal::Matrix& iwnl, const int l, const internal::Array& igrid);

         /**
          * @brief Compute \f$W_0^l (r)\f$
          *
          * Internal computation can be done in multiple precision
          */
         static void W0l(Eigen::Ref<internal::Matrix> iw0ab, const int l, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, JacobiPolynomial::NormalizerAB norm);

         /**
          * @brief Polynomial normalizer for unit Worland normalization
          */
         static internal::Array unitWPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Polynomial n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Polynomial n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief First derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWD2Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWD2P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Second derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWD2P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative normalizer for unit Worland normalization
          */
         static internal::Array unitWD3Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative n=0 normalizer for unit Worland normalization
          */
         static internal::Array unitWD3P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

         /**
          * @brief Third derivative n=1 normalizer for unit Worland normalization
          */
         static internal::Array unitWD3P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta);

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
