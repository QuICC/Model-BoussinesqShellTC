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

   JacobiPolynomial::NormalizerNAB JacobiPolynomial::mpNormPnab = &JacobiPolynomial::naturalPnab;
   JacobiPolynomial::NormalizerAB JacobiPolynomial::mpNormP0ab = &JacobiPolynomial::naturalP0ab;
   JacobiPolynomial::NormalizerAB JacobiPolynomial::mpNormP1ab = &JacobiPolynomial::naturalP1ab;

   JacobiPolynomial::NormalizerNAB JacobiPolynomial::mpNormDPnab = &JacobiPolynomial::naturalDPnab;
   JacobiPolynomial::NormalizerAB JacobiPolynomial::mpNormDP0ab = &JacobiPolynomial::naturalDP0ab;
   JacobiPolynomial::NormalizerAB JacobiPolynomial::mpNormDP1ab = &JacobiPolynomial::naturalDP1ab;

   void JacobiPolynomial::setNormalization(const JacobiPolynomial::NormalizationId id)
   {
      if(id == JacobiPolynomial::NATURAL)
      {
         JacobiPolynomial::mpNormPnab = &JacobiPolynomial::naturalPnab;
         JacobiPolynomial::mpNormP0ab = &JacobiPolynomial::naturalP0ab;
         JacobiPolynomial::mpNormP1ab = &JacobiPolynomial::naturalP1ab;

         JacobiPolynomial::mpNormDPnab = &JacobiPolynomial::naturalDPnab;
         JacobiPolynomial::mpNormDP0ab = &JacobiPolynomial::naturalDP0ab;
         JacobiPolynomial::mpNormDP1ab = &JacobiPolynomial::naturalDP1ab;

      } else if(id == JacobiPolynomial::UNITWORLAND)
      {
         JacobiPolynomial::mpNormPnab = &JacobiPolynomial::unitWPnab;
         JacobiPolynomial::mpNormP0ab = &JacobiPolynomial::unitWP0ab;
         JacobiPolynomial::mpNormP1ab = &JacobiPolynomial::unitWP1ab;

         JacobiPolynomial::mpNormDPnab = &JacobiPolynomial::unitWDPnab;
         JacobiPolynomial::mpNormDP0ab = &JacobiPolynomial::unitWDP0ab;
         JacobiPolynomial::mpNormDP1ab = &JacobiPolynomial::unitWDP1ab;
      }
   }

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
      JacobiPolynomial::P0ab(ipoly.col(0), alpha, beta, igrid, JacobiPolynomial::mpNormP0ab);

      if(nPoly > 1)
      {
         JacobiPolynomial::P1ab(ipoly.col(1), alpha, beta, ipoly.col(0), igrid, JacobiPolynomial::mpNormP1ab);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         JacobiPolynomial::Pnab(ipoly.col(i), i, alpha, beta, ipoly.col(i-1), ipoly.col(i-2), igrid, JacobiPolynomial::mpNormPnab);
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

      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         JacobiPolynomial::P0ab(idiff.col(1), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), igrid, JacobiPolynomial::mpNormDP0ab);
      }

      if(nPoly > 2)
      {
         JacobiPolynomial::P1ab(idiff.col(2), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idiff.col(1), igrid, JacobiPolynomial::mpNormDP1ab);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         JacobiPolynomial::Pnab(idiff.col(i), i-1, alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idiff.col(i-1), idiff.col(i-2), igrid, JacobiPolynomial::mpNormDPnab);
      }

      diff = Precision::cast(idiff);
   }

   void JacobiPolynomial::drrPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat e, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(alpha < -1 || beta < -1)
      {
         throw Exception("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix ipnab(gN,2);
      internal::Matrix idpnab(gN,2);

      // Compute P_0
      JacobiPolynomial::P0ab(ipnab.col(0), alpha, beta, igrid, JacobiPolynomial::mpNormP0ab);

      // Compute DP_0
      idpnab.col(0).setZero();

      // Compute e P + 2(x+1) DP
      idiff.resize(gN, nPoly);
      JacobiPolynomial::drrP0ab(idiff.col(0), e, ipnab.col(0));

      if(nPoly > 1)
      {
         // Compute P_1
         JacobiPolynomial::P1ab(ipnab.col(1), alpha, beta, ipnab.col(0), igrid, JacobiPolynomial::mpNormP1ab);

         // Compute DP_1
         JacobiPolynomial::P0ab(idpnab.col(0), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), igrid, JacobiPolynomial::mpNormDP0ab);

         // Compute e P + 2(x+1) DP
         JacobiPolynomial::drrPnab(idiff.col(1), e, ipnab.col(1), idpnab.col(0), igrid);
      }

      if(nPoly > 2)
      {
         // Increment P_n
         JacobiPolynomial::Pnab(ipnab.col(0), 2, alpha, beta, ipnab.col(1), idiff.col(0), igrid, JacobiPolynomial::mpNormPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Compute DP_2
         JacobiPolynomial::P1ab(idpnab.col(1), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idpnab.col(0), igrid, JacobiPolynomial::mpNormDP1ab);

         // Compute e P + 2(x+1) DP
         JacobiPolynomial::drrPnab(idiff.col(2), e, ipnab.col(1), idpnab.col(1), igrid);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         // Increment P_n
         JacobiPolynomial::Pnab(ipnab.col(0), i, alpha, beta, idiff.col(1), ipnab.col(0), igrid, JacobiPolynomial::mpNormPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Increment DP_n
         JacobiPolynomial::Pnab(idpnab.col(0), i-1, alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idpnab.col(1), idpnab.col(0), igrid, JacobiPolynomial::mpNormDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute e P + 2(x+1) DP
         JacobiPolynomial::drrPnab(idiff.col(i), e, ipnab.col(1), idpnab.col(1), igrid);
      }

      diff = Precision::cast(idiff);
   }

   void JacobiPolynomial::drPnab(Matrix& diff, internal::Matrix& idiff, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(alpha < -1 || beta < -1)
      {
         throw Exception("Tried to compute Jacobi polynomial P_n^{(alpha,beta)} with alpha < -1 or beta < -1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      // Storage for dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);

      // Compute DP_0
      idpnab.col(0).setZero();

      // Compute 2(x+1) DP
      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         // Compute DP_1
         JacobiPolynomial::P0ab(idpnab.col(0), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), igrid, JacobiPolynomial::mpNormDP0ab);

         // Compute 2(x+1) DP
         JacobiPolynomial::drPnab(idiff.col(1), idpnab.col(0), igrid);
      }

      if(nPoly > 2)
      {
         // Compute DP_2
         JacobiPolynomial::P1ab(idpnab.col(1), alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idpnab.col(0), igrid, JacobiPolynomial::mpNormDP1ab);

         // Compute 2(x+1) DP
         JacobiPolynomial::drPnab(idiff.col(2), idpnab.col(1), igrid);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         // Increment DP_n
         JacobiPolynomial::Pnab(idpnab.col(0), i-1, alpha + MHD_MP(1.0), beta + MHD_MP(1.0), idpnab.col(1), idpnab.col(0), igrid, JacobiPolynomial::mpNormDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute 2(x+1) DP
         JacobiPolynomial::drPnab(idiff.col(i), idpnab.col(1), igrid);
      }

      diff = Precision::cast(idiff);
   }

   void JacobiPolynomial::Pnab(Eigen::Ref<internal::Matrix> ipnab, const int n, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ipn_1ab, const Eigen::Ref<const internal::Matrix>& ipn_2ab, const internal::Array& igrid, NormalizerNAB norm)
   {
      internal::MHDFloat dn = internal::MHDFloat(n);
      internal::Array cs = norm(dn, alpha, beta);

      ipnab.array() = cs(0)*ipn_2ab.array();
      ipnab.array() += (cs(1)*igrid.array() + cs(2))*ipn_1ab.array();
      ipnab.array() *= cs(3);
   }

   void JacobiPolynomial::P1ab(Eigen::Ref<internal::Matrix> ip1ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const Eigen::Ref<const internal::Matrix>& ip0ab, const internal::Array& igrid, NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      ip1ab.array() = (cs(0)*igrid.array() + cs(1))*ip0ab.array();
      ip1ab.array() *= cs(2);
   }

   void JacobiPolynomial::P0ab(Eigen::Ref<internal::Matrix> ip0ab, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      ip0ab.setConstant(cs(0));
   }

   void JacobiPolynomial::drrPnab(Eigen::Ref<internal::Matrix> idrrpnab, const internal::MHDFloat e, const Eigen::Ref<const internal::Matrix>& ipnab, const Eigen::Ref<const internal::Matrix>& idpnab, const internal::Array& igrid)
   {
      idrrpnab.array() = e*ipnab.array() + MHD_MP(2.0)*(igrid.array() + MHD_MP(1.0))*idpnab.array();
   }

   void JacobiPolynomial::drrP0ab(Eigen::Ref<internal::Matrix> idrrp0ab, const internal::MHDFloat e, const Eigen::Ref<const internal::Matrix>& ipnab)
   {
      idrrp0ab.array() = e*ipnab;
   }

   void JacobiPolynomial::drPnab(Eigen::Ref<internal::Matrix> idrpnab, const Eigen::Ref<const internal::Matrix>& idpnab, const internal::Array& igrid)
   {
      idrpnab.array() = MHD_MP(2.0)*(igrid.array() + MHD_MP(1.0))*idpnab.array();
   }

   //
   // Natural polynomial normalizer
   //
   internal::Array JacobiPolynomial::naturalPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta)/MHD_MP(2.0);
      cs(1) = (alpha - beta)/MHD_MP(2.0);
      cs(2) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   //
   // Natural polynomial normalizer
   //
   internal::Array JacobiPolynomial::naturalDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta)/MHD_MP(2.0);
      cs(1) = (alpha - beta)/MHD_MP(2.0);
      cs(2) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::naturalDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   //
   // Unit Worland polynomial normalizers
   //

   internal::Array JacobiPolynomial::unitWPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::unitWP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta)/MHD_MP(2.0);
      cs(1) = (alpha - beta)/MHD_MP(2.0);
      cs(2) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::unitWP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

      return cs;
   }

   //
   // Unit Worland derivative normalizers
   //
   internal::Array JacobiPolynomial::unitWDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::unitWDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta)/MHD_MP(2.0);
      cs(1) = (alpha - beta)/MHD_MP(2.0);
      cs(2) = MHD_MP(1.0);

      return cs;
   }

   internal::Array JacobiPolynomial::unitWDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(1.0);

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
