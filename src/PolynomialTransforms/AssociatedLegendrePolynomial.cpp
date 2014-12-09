/** 
 * @file AssociatedLegendrePolynomial.cpp
 * @brief Source of the implementation of the associated Legendre polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/AssociatedLegendrePolynomial.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Quadratures/LegendreRule.hpp"

namespace GeoMHDiSCC {

namespace Polynomial {

   void AssociatedLegendrePolynomial::Plm(Matrix& op, const int m, const internal::Array& igrid)
   {
   }

   void AssociatedLegendrePolynomial::dPlm(Matrix& op, const int m, const internal::Array& igrid)
   {
   }

   void AssociatedLegendrePolynomial::Pmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op.setConstant(MHD_MP(1.0));
      } else
      {
         internal::MHDFloat di;
         internal::MHDFloat factor = MHD_MP(1.0);

         for(int i = 1; i <= m; i++)
         {
            di = internal::MHDFloat(2*i);
            factor *= -precision::sqrt((di - MHD_MP(1.0))/di);
         }
         op = igrid.array().acos();
         op = op.array().sin();
         op = op.array().pow(m);
         op *= factor;
      }
   }

   void AssociatedLegendrePolynomial::Pmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op *= precision::sqrt(internal::MHDFloat(2*m + 1));
         op.array() *= ipmm.array(); 
      }
   }

   void AssociatedLegendrePolynomial::dPmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
   }

   void AssociatedLegendrePolynomial::dPmm1(internal::Array& op, const int m, const internal::Array& igrid)
   {
   }

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial()
   {
   }

   AssociatedLegendrePolynomial::~AssociatedLegendrePolynomial()
   {
   }

}
}
