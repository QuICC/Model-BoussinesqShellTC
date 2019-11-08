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

namespace QuICC {

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
      ThreeTermRecurrence::P0(ipoly.col(0), alpha, beta, igrid, JacobiPolynomial::normP0ab());

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), alpha, beta, ipoly.col(0), igrid, JacobiPolynomial::normP1ab());
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(i), i, alpha, beta, ipoly.col(i-1), ipoly.col(i-2), igrid, JacobiPolynomial::normPnab());
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
         ThreeTermRecurrence::P0(idiff.col(1), a1, b1, igrid, JacobiPolynomial::normDP0ab());
      }

      if(nPoly > 2)
      {
         ThreeTermRecurrence::P1(idiff.col(2), a1, b1, idiff.col(1), igrid, JacobiPolynomial::normDP1ab());
      }

      for(int i = 3; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-1, a1, b1, idiff.col(i-1), idiff.col(i-2), igrid, JacobiPolynomial::normDPnab());
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
         ThreeTermRecurrence::P0(idiff.col(2), a2, b2, igrid, JacobiPolynomial::normD2P0ab());
      }

      if(nPoly > 3)
      {
         ThreeTermRecurrence::P1(idiff.col(3), a2, b2, idiff.col(2), igrid, JacobiPolynomial::normD2P1ab());
      }

      for(int i = 4; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-2, a2, b2, idiff.col(i-1), idiff.col(i-2), igrid, JacobiPolynomial::normD2Pnab());
      }

      diff = Precision::cast(idiff);
   }

   void JacobiPolynomial::d3Pnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
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

      internal::MHDFloat a3 = alpha + MHD_MP(3.0);
      internal::MHDFloat b3 = beta + MHD_MP(3.0);

      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         idiff.col(1).setZero();
      }

      if(nPoly > 2)
      {
         idiff.col(2).setZero();
      }

      if(nPoly > 3)
      {
         ThreeTermRecurrence::P0(idiff.col(3), a3, b3, igrid, JacobiPolynomial::normD3P0ab());
      }

      if(nPoly > 4)
      {
         ThreeTermRecurrence::P1(idiff.col(4), a3, b3, idiff.col(3), igrid, JacobiPolynomial::normD3P1ab());
      }

      for(int i = 5; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-3, a3, b3, idiff.col(i-1), idiff.col(i-2), igrid, JacobiPolynomial::normD3Pnab());
      }

      diff = Precision::cast(idiff);
   }

   //
   // General polynomial normalizer
   //
   ThreeTermRecurrence::NormalizerNAB  JacobiPolynomial::normPnab()
   {
      return &JacobiPolynomial::naturalPnab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normP1ab()
   {
      return &JacobiPolynomial::naturalP1ab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normP0ab()
   {
      return &JacobiPolynomial::naturalP0ab;
   }

   ThreeTermRecurrence::NormalizerNAB  JacobiPolynomial::normDPnab()
   {
      return &JacobiPolynomial::naturalDPnab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normDP1ab()
   {
      return &JacobiPolynomial::naturalDP1ab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normDP0ab()
   {
      return &JacobiPolynomial::naturalDP0ab;
   }

   ThreeTermRecurrence::NormalizerNAB  JacobiPolynomial::normD2Pnab()
   {
      return &JacobiPolynomial::naturalD2Pnab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normD2P1ab()
   {
      return &JacobiPolynomial::naturalD2P1ab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normD2P0ab()
   {
      return &JacobiPolynomial::naturalD2P0ab;
   }

   ThreeTermRecurrence::NormalizerNAB  JacobiPolynomial::normD3Pnab()
   {
      return &JacobiPolynomial::naturalD3Pnab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normD3P1ab()
   {
      return &JacobiPolynomial::naturalD3P1ab;
   }

   ThreeTermRecurrence::NormalizerAB  JacobiPolynomial::normD3P0ab()
   {
      return &JacobiPolynomial::naturalD3P0ab;
   }

   //
   // Natural polynomial normalizer
   //
   internal::Array JacobiPolynomial::naturalPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n*(n + a + b));
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(n + a + b)*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = MHD_MP(0.5);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural first derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalDPnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(2.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(1.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.5)*(a + b);

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural second derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalD2Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(3.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(2.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalD2P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(1.0)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalD2P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.25)*(a + b)*(a + b - MHD_MP(1.0));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   //
   // Natural third derivative normalizer
   //
   internal::Array JacobiPolynomial::naturalD3Pnab(const internal::MHDFloat n, const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(4);

      cs(0) = -((n + a + b - MHD_MP(1.0))*(n + a - MHD_MP(1.0))*(n + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(n*(n + a + b - MHD_MP(4.0))*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(MHD_MP(2.0)*n + a + b))/(MHD_MP(2.0)*n);
      cs(2) = ((MHD_MP(2.0)*n + a + b - MHD_MP(1.0))*(a*a - b*b))/(MHD_MP(2.0)*n*(MHD_MP(2.0)*n + a + b - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(n + a + b - MHD_MP(3.0));

      assert(!precision::isnan(cs.sum()));
      
      return cs;
   }

   internal::Array JacobiPolynomial::naturalD3P1ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(3);

      cs(0) = (a + b + MHD_MP(2.0));
      cs(1) = (a - b);
      cs(2) = (a + b + MHD_MP(1.0))/(MHD_MP(2.0)*(a + b - MHD_MP(2.0)));

      assert(!precision::isnan(cs.sum()));

      return cs;
   }

   internal::Array JacobiPolynomial::naturalD3P0ab(const internal::MHDFloat a, const internal::MHDFloat b)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.125)*(a + b)*(a + b - MHD_MP(1.0))*(a + b - MHD_MP(2.0));

      assert(!precision::isnan(cs.sum()));

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
