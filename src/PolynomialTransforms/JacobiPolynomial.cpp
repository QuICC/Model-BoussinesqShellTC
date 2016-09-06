/** 
 * @file JacobiPolynomial.cpp
 * @brief Source of the implementation of the Jacobi polynomial
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "PolynomialTransforms/JacobiPolynomial.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace Polynomial {

   void JacobiPolynomial::Pnab(Matrix& poly, internal::Matrix& ipoly, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (alpha < -1 || beta < -1)
      {
         throw Exception("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      ipoly.resize(gN, nPoly);
      ThreeTermRecurrence::P0(ipoly.col(0), alpha, beta, igrid, &JacobiPolynomial::naturalP0ab);

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), alpha, beta, ipoly.col(0), igrid, &JacobiPolynomial::naturalP1ab);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(i), i, alpha, beta, ipoly.col(i-1), ipoly.col(i-2), igrid, &JacobiPolynomial::naturalPnab);
      }

      poly = Precision::cast(ipoly);
   }

   void JacobiPolynomial::dPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (alpha < -1 || beta < -1)
      {
         throw Exception("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      internal::MHDFloat a1 = alpha + MHD_MP(1.0);
      internal::MHDFloat b1 = beta + MHD_MP(1.0);

      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P0(idiff.col(1), a1, b1, igrid, &JacobiPolynomial::naturalDP0ab);
      }

      if(nPoly > 2)
      {
         ThreeTermRecurrence::P1(idiff.col(2), a1, b1, idiff.col(1), igrid, &JacobiPolynomial::naturalDP1ab);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-1, a1, b1, idiff.col(i-1), idiff.col(i-2), igrid, &JacobiPolynomial::naturalDPnab);
      }

      diff = Precision::cast(idiff);
   }

   void JacobiPolynomial::d2Pnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (alpha < -1 || beta < -1)
      {
         throw Exception("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      internal::MHDFloat a2 = alpha + MHD_MP(2.0);
      internal::MHDFloat b2 = beta + MHD_MP(2.0);

      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         idiff.col(1).setZero();
      }

      if(nPoly > 2)
      {
         ThreeTermRecurrence::P0(idiff.col(2), a2, b2, igrid, &JacobiPolynomial::naturalD2P0ab);
      }

      if(nPoly > 3)
      {
         ThreeTermRecurrence::P1(idiff.col(3), a2, b2, idiff.col(2), igrid, &JacobiPolynomial::naturalD2P1ab);
      }

      for(int i = 4; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-2, a2, b2, idiff.col(i-1), idiff.col(i-2), igrid, &JacobiPolynomial::naturalD2Pnab);
      }

      diff = Precision::cast(idiff);
   }

   //
   // Natural polynomial normalizer
   //
   internal::Array JacobiPolynomial::naturalPnab(const internal::MHDFloat n, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((n + alpha - MHD_MP(1.0))*(n + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(n*(n + alpha + beta)*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(MHD_MP(2.0)*n*(n + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*n*(n + alpha + beta)*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (alpha + beta + MHD_MP(2.0));
      cs(1) = (alpha - beta);
      cs(2) = MHD_MP(0.5);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   //
   // Natural first derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalDPnab(const internal::MHDFloat n, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((n + alpha + beta - MHD_MP(1.0))/(n + alpha + beta - MHD_MP(2.0)))*((n + alpha - MHD_MP(1.0))*(n + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + alpha + beta - MHD_MP(1.0));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (alpha + beta + MHD_MP(2.0));
      cs(1) = (alpha - beta);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*(alpha + beta));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = (alpha + beta);

      return cs;
   }

   //
   // Natural second derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalD2Pnab(const internal::MHDFloat n, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((n + alpha + beta - MHD_MP(1.0))*(n + alpha - MHD_MP(1.0))*(n + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/((n + alpha + beta - MHD_MP(3.0))*n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + alpha + beta - MHD_MP(2.0));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalD2P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta);
      cs(1) = (alpha - beta);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*(alpha + beta - MHD_MP(1.0)));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalD2P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.25)*precision::exp(precisiontr1::lgamma(alpha + beta + MHD_MP(1.0)) - precisiontr1::lgamma(alpha + beta - MHD_MP(1.0)));

      return cs;
   }

   //
   // Natural third derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalD3Pnab(const internal::MHDFloat n, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((n + alpha + beta - MHD_MP(1.0))*(n + alpha - MHD_MP(1.0))*(n + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/((n + alpha + beta - MHD_MP(4.0))*n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*n + alpha + beta))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + alpha + beta - MHD_MP(3.0));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalD3P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta);
      cs(1) = (alpha - beta);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*(alpha + beta - MHD_MP(1.0)));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalD3P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.125)*precision::exp(precisiontr1::lgamma(alpha + beta + MHD_MP(1.0)) - precisiontr1::lgamma(alpha + beta - MHD_MP(2.0)));

      return cs;
   }

   JacobiPolynomial::JacobiPolynomial()
   {
   }

   JacobiPolynomial::~JacobiPolynomial()
   {
   }

}
}
