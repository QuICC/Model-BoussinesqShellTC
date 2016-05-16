/** 
 * @file WorlandPolynomial.cpp
 * @brief Source of the implementation of the Jones-Worland polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/WorlandPolynomial.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "PolynomialTransforms/JacobiPolynomial.hpp"

namespace GeoMHDiSCC {

namespace Polynomial {

   internal::MHDFloat WorlandPolynomial::alpha(const int l)
   {
      return MHD_MP(-0.5);
   }

   internal::MHDFloat WorlandPolynomial::beta(const int l)
   {
      return internal::MHDFloat(l) - MHD_MP(0.5);
   }

   void WorlandPolynomial::Wnl(Matrix& poly, internal::Matrix& ipoly, const int l, const internal::Array& igrid)
   {
      if (l < 0)
      {
         throw Exception("Tried to compute Worland polynomial W_n^l with l < 0");
      }

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      JacobiPolynomial::Pnab(poly, ipoly, WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ixgrid);

      WorlandPolynomial::rl(ipoly, l, igrid);

      poly = Precision::cast(ipoly);
   }

   void WorlandPolynomial::dWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      if (l < 0)
      {
         throw Exception("Tried to compute Worland polynomial W_n^l with l < 0");
      }

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      if(l > 0)
      {
         JacobiPolynomial::drrPnab(diff, idiff, internal::MHDFloat(l), WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ixgrid);

         WorlandPolynomial::rl(idiff, l-1, igrid);
      } else
      {
         JacobiPolynomial::drPnab(diff, idiff, WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ixgrid);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::Wnl_r(Matrix& poly, internal::Matrix& ipoly, const int l, const internal::Array& igrid)
   {
      if (l < 0)
      {
         throw Exception("Tried to compute Worland polynomial W_n^l with l < 0");
      }

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      JacobiPolynomial::Pnab(poly, ipoly, WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ixgrid);

      WorlandPolynomial::rl(ipoly, l-1, igrid);

      poly = Precision::cast(ipoly);
   }

   void WorlandPolynomial::rl(internal::Matrix& ipoly, const int l, const internal::Array& igrid)
   {
      if (l < 0)
      {
         throw Exception("Tried to multiply polynomial with 1/r^a with a > 0!");
      }

      if(l > 0)
      {
         internal::Array factor = igrid.array().pow(l);

         for(int i = 0; i < ipoly.cols(); i++)
         {
            ipoly.col(i).array() *= factor.array();
         }
      }
   }

   WorlandPolynomial::WorlandPolynomial()
   {
   }

   WorlandPolynomial::~WorlandPolynomial()
   {
   }

}
}
