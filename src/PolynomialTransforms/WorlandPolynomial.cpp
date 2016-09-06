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
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (l < 0)
      {
         throw Exception("Tried to compute Worland polynomial W_n^l with l < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a = WorlandPolynomial::alpha(l);
      internal::MHDFloat b = WorlandPolynomial::beta(l);

      ipoly.resize(gN, nPoly);
      WorlandPolynomial::W0l(ipoly.col(0), l, a, b, igrid, &WorlandPolynomial::unitWP0ab);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), a, b, ipoly.col(0), ixgrid, &WorlandPolynomial::unitWP1ab);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(i), i, a, b, ipoly.col(i-1), ipoly.col(i-2), ixgrid, &WorlandPolynomial::unitWPnab);
      }

      poly = Precision::cast(ipoly);
   }

   void WorlandPolynomial::dWn0(Matrix& diff, internal::Matrix& idiff, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = WorlandPolynomial::alpha(0) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(0) + MHD_MP(1.0);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      idiff.resize(gN, nPoly);
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         WorlandPolynomial::W0l(idiff.col(1), 1, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idiff.col(1) *= MHD_MP(4.0);
      }

      if(nPoly > 2)
      {
         ThreeTermRecurrence::P1(idiff.col(2), a1, b1, idiff.col(1), ixgrid, &WorlandPolynomial::unitWDP1ab);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(idiff.col(i), i-1, a1, b1, idiff.col(i-1), idiff.col(i-2), ixgrid, &WorlandPolynomial::unitWDPnab);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::dWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      if(l == 0)
      {
         WorlandPolynomial::dWn0(diff, idiff, igrid);
      } else
      {
         int gN = diff.rows();
         int nPoly = diff.cols();

         if (nPoly < 1)
         {
            throw Exception("Operator matrix should have at least 1 column");
         }

         if (gN != igrid.size())
         {
            throw Exception("Operator matrix does not mach grid size");
         }

         internal::MHDFloat a = WorlandPolynomial::alpha(l);
         internal::MHDFloat b = WorlandPolynomial::beta(l);
         internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
         internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
         internal::MHDFloat dl = internal::MHDFloat(l);

         // Make X grid in [-1, 1]
         internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

         // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
         internal::Matrix ipnab(gN,2);
         internal::Matrix idpnab(gN,2);
         idiff.resize(gN, nPoly);

         // Compute P_0
         WorlandPolynomial::W0l(ipnab.col(0), l-1, a, b, igrid, &WorlandPolynomial::unitWP0ab);
         ipnab.col(0) *= dl;

         // Compute DP_0
         idpnab.col(0).setZero();

         // Compute l P
         idiff.col(0) = ipnab.col(0);

         if(nPoly > 1)
         {
            // Compute P_0
            ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, &WorlandPolynomial::unitWP1ab);

            // Compute DP_1
            WorlandPolynomial::W0l(idpnab.col(0), l+1, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
            idpnab.col(0) *= MHD_MP(4.0);

            // Compute e P + 4r^2 DP
            idiff.col(1) = ipnab.col(1) + idpnab.col(0);
         }

         if(nPoly > 2)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
            ipnab.col(0).swap(ipnab.col(1));

            // Compute DP_2
            ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

            // Compute e P + 2(x+1) DP
            idiff.col(2) = ipnab.col(1) + idpnab.col(1);
         }

         for(int i = 3; i < nPoly; ++i)
         {
            // Increment P_n
            ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
            ipnab.col(0).swap(ipnab.col(1));

            // Increment DP_n
            ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
            idpnab.col(0).swap(idpnab.col(1));

            // Compute e P + 2(x+1) DP
            idiff.col(i) = ipnab.col(1) + idpnab.col(1);
         }

         diff = Precision::cast(idiff);
      }
   }

   void WorlandPolynomial::drWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a = WorlandPolynomial::alpha(l);
      internal::MHDFloat b = WorlandPolynomial::beta(l);
      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat dl1 = internal::MHDFloat(l+1);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix ipnab(gN,2);
      internal::Matrix idpnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute P_0
      WorlandPolynomial::W0l(ipnab.col(0), l, a, b, igrid, &WorlandPolynomial::unitWP0ab);
      ipnab.col(0) *= dl1;

      // Compute DP_0
      idpnab.col(0).setZero();

      // Compute l P
      idiff.col(0) = ipnab.col(0);

      if(nPoly > 1)
      {
         // Compute P_0
         ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, &WorlandPolynomial::unitWP1ab);

         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l+2, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(4.0);

         // Compute e P + 4r^2 DP
         idiff.col(1) = ipnab.col(1) + idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment P_n
         ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Compute DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(2) = ipnab.col(1) + idpnab.col(1);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         // Increment P_n
         ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = ipnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::r_1drWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(l < 1)
      {
         throw Exception("Tried to compute Worland polynomial 1/r d/dr r W_n^l with l < 1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a = WorlandPolynomial::alpha(l);
      internal::MHDFloat b = WorlandPolynomial::beta(l);
      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat dl1 = internal::MHDFloat(l+1);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix ipnab(gN,2);
      internal::Matrix idpnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute P_0
      WorlandPolynomial::W0l(ipnab.col(0), l-1, a, b, igrid, &WorlandPolynomial::unitWP0ab);
      ipnab.col(0) *= dl1;

      // Compute DP_0
      idpnab.col(0).setZero();

      // Compute l P
      idiff.col(0) = ipnab.col(0);

      if(nPoly > 1)
      {
         // Compute P_0
         ThreeTermRecurrence::P1(ipnab.col(1), a, b, ipnab.col(0), ixgrid, &WorlandPolynomial::unitWP1ab);

         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l+1, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(4.0);

         // Compute e P + 4r^2 DP
         idiff.col(1) = ipnab.col(1) + idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment P_n
         ThreeTermRecurrence::Pn(ipnab.col(0), 2, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Compute DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(2) = ipnab.col(1) + idpnab.col(1);
      }

      for(int i = 3; i < nPoly; ++i)
      {
         // Increment P_n
         ThreeTermRecurrence::Pn(ipnab.col(0), i, a, b, ipnab.col(1), ipnab.col(0), ixgrid, &WorlandPolynomial::unitWPnab);
         ipnab.col(0).swap(ipnab.col(1));

         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = ipnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::slaplWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(l < 0)
      {
         throw Exception("Tried to compute Worland spherical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat a2 = WorlandPolynomial::alpha(l) + MHD_MP(2.0);
      internal::MHDFloat b2 = WorlandPolynomial::beta(l) + MHD_MP(2.0);
      internal::MHDFloat dl = internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);
      internal::Matrix id2pnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute spherical laplacian P_0
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(4.0)*(MHD_MP(2.0)*dl + MHD_MP(3.0)); 

         // Compute spherical laplacian P_1
         idiff.col(1) = idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute D2P_2
         WorlandPolynomial::W0l(id2pnab.col(0), l+2, a2, b2, igrid, &WorlandPolynomial::unitWD2P0ab);
         id2pnab.col(0) *= MHD_MP(16.0); 

         // Compute e P + 2(x+1) DP
         idiff.col(2) = id2pnab.col(0) + idpnab.col(1);
      }

      if(nPoly > 3)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2P1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(3) = id2pnab.col(1) + idpnab.col(1);
      }

      for(int i = 4; i < nPoly; ++i)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2Pnab);
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = id2pnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::claplhWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(l < 0)
      {
         throw Exception("Tried to compute Worland spherical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat a2 = WorlandPolynomial::alpha(l) + MHD_MP(2.0);
      internal::MHDFloat b2 = WorlandPolynomial::beta(l) + MHD_MP(2.0);
      internal::MHDFloat dl = internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);
      internal::Matrix id2pnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute spherical laplacian P_0
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(8.0)*(dl + MHD_MP(1.0)); 

         // Compute spherical laplacian P_1
         idiff.col(1) = idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute D2P_2
         WorlandPolynomial::W0l(id2pnab.col(0), l+2, a2, b2, igrid, &WorlandPolynomial::unitWD2P0ab);
         id2pnab.col(0) *= MHD_MP(16.0); 

         // Compute e P + 2(x+1) DP
         idiff.col(2) = id2pnab.col(0) + idpnab.col(1);
      }

      if(nPoly > 3)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2P1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(3) = id2pnab.col(1) + idpnab.col(1);
      }

      for(int i = 4; i < nPoly; ++i)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2Pnab);
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = id2pnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::r_1claplhWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(l < 0)
      {
         throw Exception("Tried to compute Worland spherical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat a2 = WorlandPolynomial::alpha(l) + MHD_MP(2.0);
      internal::MHDFloat b2 = WorlandPolynomial::beta(l) + MHD_MP(2.0);
      internal::MHDFloat dl = internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);
      internal::Matrix id2pnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute spherical laplacian P_0
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l-1, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(8.0)*(dl + MHD_MP(1.0)); 

         // Compute spherical laplacian P_1
         idiff.col(1) = idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute D2P_2
         WorlandPolynomial::W0l(id2pnab.col(0), l+1, a2, b2, igrid, &WorlandPolynomial::unitWD2P0ab);
         id2pnab.col(0) *= MHD_MP(16.0); 

         // Compute e P + 2(x+1) DP
         idiff.col(2) = id2pnab.col(0) + idpnab.col(1);
      }

      if(nPoly > 3)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2P1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(3) = id2pnab.col(1) + idpnab.col(1);
      }

      for(int i = 4; i < nPoly; ++i)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2Pnab);
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = id2pnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::dclaplhWnl(Matrix& diff, internal::Matrix& idiff, const int l, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if(l < 0)
      {
         throw Exception("Tried to compute Worland spherical laplacian with l < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      internal::MHDFloat a1 = WorlandPolynomial::alpha(l) + MHD_MP(1.0);
      internal::MHDFloat b1 = WorlandPolynomial::beta(l) + MHD_MP(1.0);
      internal::MHDFloat a2 = WorlandPolynomial::alpha(l) + MHD_MP(2.0);
      internal::MHDFloat b2 = WorlandPolynomial::beta(l) + MHD_MP(2.0);
      internal::MHDFloat a3 = WorlandPolynomial::alpha(l) + MHD_MP(3.0);
      internal::MHDFloat b3 = WorlandPolynomial::beta(l) + MHD_MP(3.0);
      internal::MHDFloat dl = internal::MHDFloat(l);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      // Storage for P_n^{(alpha,beta)} and dP_n{(alpha,beta)}
      internal::Matrix idpnab(gN,2);
      internal::Matrix id2pnab(gN,2);
      internal::Matrix id3pnab(gN,2);
      idiff.resize(gN, nPoly);

      // Compute spherical laplacian P_0
      idiff.col(0).setZero();

      if(nPoly > 1)
      {
         // Compute DP_1
         WorlandPolynomial::W0l(idpnab.col(0), l-1, a1, b1, igrid, &WorlandPolynomial::unitWDP0ab);
         idpnab.col(0) *= MHD_MP(8.0)*dl*(dl + MHD_MP(1.0)); 

         // Compute spherical laplacian P_1
         idiff.col(1) = idpnab.col(0);
      }

      if(nPoly > 2)
      {
         // Increment DP_2
         ThreeTermRecurrence::P1(idpnab.col(1), a1, b1, idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDP1ab);

         // Compute D2P_2
         WorlandPolynomial::W0l(id2pnab.col(0), l+1, a2, b2, igrid, &WorlandPolynomial::unitWD2P0ab);
         id2pnab.col(0) *= MHD_MP(16.0)*(MHD_MP(3.0)*dl + MHD_MP(4.0)); 

         // Compute e P + 2(x+1) DP
         idiff.col(2) = id2pnab.col(0) + idpnab.col(1);
      }

      if(nPoly > 3)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 2, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::P1(id2pnab.col(1), a2, b2, id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2P1ab);

         // Compute D3P_2
         WorlandPolynomial::W0l(id3pnab.col(0), l+3, a3, b3, igrid, &WorlandPolynomial::unitWD3P0ab);
         id2pnab.col(0) *= MHD_MP(64.0); 

         // Compute e P + 2(x+1) DP
         idiff.col(3) = id2pnab.col(0) + id2pnab.col(1) + idpnab.col(1);
      }

      if(nPoly > 4)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), 3, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Compute D2P_2
         ThreeTermRecurrence::Pn(id2pnab.col(0), 2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2Pnab);
         id2pnab.col(0).swap(id2pnab.col(1));

         // Compute D3P_2
         ThreeTermRecurrence::P1(id3pnab.col(1), a3, b3, id3pnab.col(0), ixgrid, &WorlandPolynomial::unitWD3P1ab);

         // Compute e P + 2(x+1) DP
         idiff.col(4) = id3pnab.col(1) + id2pnab.col(1) + idpnab.col(1);
      }

      for(int i = 5; i < nPoly; ++i)
      {
         // Increment DP_n
         ThreeTermRecurrence::Pn(idpnab.col(0), i-1, a1, b1, idpnab.col(1), idpnab.col(0), ixgrid, &WorlandPolynomial::unitWDPnab);
         idpnab.col(0).swap(idpnab.col(1));

         // Increment D2P_n
         ThreeTermRecurrence::Pn(id2pnab.col(0), i-2, a2, b2, id2pnab.col(1), id2pnab.col(0), ixgrid, &WorlandPolynomial::unitWD2Pnab);
         id2pnab.col(0).swap(id2pnab.col(1));

         // Increment D3P_n
         ThreeTermRecurrence::Pn(id3pnab.col(0), i-3, a3, b3, id3pnab.col(1), id3pnab.col(0), ixgrid, &WorlandPolynomial::unitWD3Pnab);
         id3pnab.col(0).swap(id3pnab.col(1));

         // Compute e P + 2(x+1) DP
         idiff.col(i) = id3pnab.col(1) + id2pnab.col(1) + idpnab.col(1);
      }

      diff = Precision::cast(idiff);
   }

   void WorlandPolynomial::r_1Wnl(Matrix& poly, internal::Matrix& ipoly, const int l, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (l < 1)
      {
         throw Exception("Tried to compute Worland polynomial W_n^l/r with l < 1");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      if (gN != igrid.size())
      {
         throw Exception("Operator matrix does not mach grid size");
      }

      ipoly.resize(gN, nPoly);
      WorlandPolynomial::W0l(ipoly.col(0), l-1, WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), igrid, &WorlandPolynomial::unitWP0ab);

      // Make X grid in [-1, 1]
      internal::Array ixgrid = MHD_MP(2.0)*igrid.array()*igrid.array() - MHD_MP(1.0);

      if(nPoly > 1)
      {
         ThreeTermRecurrence::P1(ipoly.col(1), WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ipoly.col(0), ixgrid, &WorlandPolynomial::unitWP1ab);
      }

      for(int i = 2; i < nPoly; ++i)
      {
         ThreeTermRecurrence::Pn(ipoly.col(i), i, WorlandPolynomial::alpha(l), WorlandPolynomial::beta(l), ipoly.col(i-1), ipoly.col(i-2), ixgrid, &WorlandPolynomial::unitWPnab);
      }

      poly = Precision::cast(ipoly);
   }

   void WorlandPolynomial::W0l(Eigen::Ref<internal::Matrix> iw0l, const int l, const internal::MHDFloat alpha, const internal::MHDFloat beta, const internal::Array& igrid, ThreeTermRecurrence::NormalizerAB norm)
   {
      internal::Array cs = norm(alpha, beta);

      if(l > 0)
      {
         iw0l.array() = igrid.array().pow(l);
      } else
      {
         iw0l.setConstant(MHD_MP(1.0));
      }

      iw0l.array() *= cs(0);
   }

   //
   // Unit Worland polynomial normalizers
   //

   internal::Array WorlandPolynomial::unitWPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta));
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(dn + alpha + beta)*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0);

      if(dn == MHD_MP(2.0))
      {
         cs(0) *= precision::sqrt((2.0*dn+alpha+beta+MHD_MP(1.0)))*precision::sqrt(((dn + alpha + beta)*dn*(dn - MHD_MP(1.0)))/((dn + alpha)*(dn + beta)*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))));
      } else
      {
         cs(0) *= precision::sqrt((2.0*dn+alpha+beta+MHD_MP(1.0))/(2.0*dn+alpha+beta-MHD_MP(3.0)))*precision::sqrt(((dn + alpha + beta)*dn*(dn + alpha + beta - MHD_MP(1.0))*(dn - MHD_MP(1.0)))/((dn + alpha)*(dn + beta)*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))));
      }
      cs(1) *= precision::sqrt((2.0*dn+alpha+beta+MHD_MP(1.0))/(2.0*dn+alpha+beta-MHD_MP(1.0)))*precision::sqrt(((dn + alpha + beta)*dn)/((dn + alpha)*(dn + beta)));
      cs(2) *= precision::sqrt((2.0*dn+alpha+beta+MHD_MP(1.0))/(2.0*dn+alpha+beta-MHD_MP(1.0)))*precision::sqrt(((dn + alpha + beta)*dn)/((dn + alpha)*(dn + beta)));
      return cs;
   }

   internal::Array WorlandPolynomial::unitWP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta);
      cs(1) = (alpha - beta);
      cs(2) = precision::sqrt((alpha + beta + MHD_MP(3.0))/(MHD_MP(4.0)*(alpha + MHD_MP(1.0))*(beta + MHD_MP(1.0))));

      return cs;
   }

   internal::Array WorlandPolynomial::unitWP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = precision::sqrt(MHD_MP(2.0))*precision::exp(MHD_MP(0.5)*(precisiontr1::lgamma(alpha + beta + MHD_MP(2.0)) - precisiontr1::lgamma(alpha + MHD_MP(1.0)) - precisiontr1::lgamma(beta + MHD_MP(1.0))));

      return cs;
   }

   //
   // Unit Worland first derivative normalizers
   //
   internal::Array WorlandPolynomial::unitWDPnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((alpha + beta + dn - MHD_MP(1.0))/(alpha + beta + dn - MHD_MP(2.0)))*((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn);
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(alpha + beta + dn - MHD_MP(1.0));

      cs(0) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(3.0)))*precision::sqrt(dn*(dn + alpha + beta - MHD_MP(2.0))/((dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))));
      cs(1) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(2) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(3) *= precision::sqrt((dn + MHD_MP(1.0))*(dn + alpha + beta - MHD_MP(1.0))/((dn + alpha)*(dn + beta)));

      return cs;
   }

   internal::Array WorlandPolynomial::unitWDP1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta)/MHD_MP(2.0);
      cs(1) = (alpha - beta)/MHD_MP(2.0);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(alpha + beta);

      cs(2) *= precision::sqrt((alpha + beta + MHD_MP(3.0))/(alpha + beta + MHD_MP(1.0)))*precision::sqrt(MHD_MP(2.0)*(alpha+beta)/((alpha + MHD_MP(1.0))*(beta + MHD_MP(1.0))));

      return cs;
   }

   internal::Array WorlandPolynomial::unitWDP0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.5)*precision::exp(precisiontr1::lgamma(alpha + beta + MHD_MP(1.0)) - precisiontr1::lgamma(alpha + beta));

      cs(0) *= precision::sqrt(alpha + beta + MHD_MP(1.0))*precision::exp(MHD_MP(0.5)*(precisiontr1::lgamma(alpha + beta) - precisiontr1::lgamma(alpha + MHD_MP(1.0)) - precisiontr1::lgamma(beta + MHD_MP(1.0))));

      return cs;
   }

   //
   // Unit Worland second derivative normalizers
   //
   internal::Array WorlandPolynomial::unitWD2Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha + beta - MHD_MP(1.0))*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/((dn + alpha + beta - MHD_MP(3.0))*dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn);
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(dn + alpha + beta - MHD_MP(2.0));

      cs(0) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))*(dn + alpha + beta - MHD_MP(3.0))*(dn + MHD_MP(1.0))/((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(3.0))*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))));
      cs(1) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(2) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(3) *= precision::sqrt((dn + alpha + beta - MHD_MP(2.0))*(dn + MHD_MP(2.0))/((dn + alpha)*(dn + beta)));

      return cs;
   }

   internal::Array WorlandPolynomial::unitWD2P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta);
      cs(1) = (alpha - beta);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*(alpha + beta - MHD_MP(1.0)));

      cs(2) *= precision::sqrt((alpha + beta + MHD_MP(3.0))*(alpha + beta - MHD_MP(1.0))*(MHD_MP(3.0))/((alpha + beta + MHD_MP(1.0))*(alpha + MHD_MP(1.0))*(beta + MHD_MP(1.0))));
      return cs;
   }

   internal::Array WorlandPolynomial::unitWD2P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.25)*precision::exp(precisiontr1::lgamma(alpha + beta + MHD_MP(1.0)) - precisiontr1::lgamma(alpha + beta - MHD_MP(1.0)));

      cs(0) *= precision::sqrt(alpha + beta + MHD_MP(1.0))*precision::exp(MHD_MP(0.5)*(precisiontr1::lgamma(alpha + beta - MHD_MP(1.0)) + precisiontr1::lgamma(MHD_MP(3.0)) - precisiontr1::lgamma(alpha + MHD_MP(1.0)) - precisiontr1::lgamma(beta + MHD_MP(1.0))));

      return cs;
   }

   //
   // Unit Worland third derivative normalizers
   //
   internal::Array WorlandPolynomial::unitWD3Pnab(const internal::MHDFloat dn, const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(4);

      cs(0) = -((dn + alpha + beta - MHD_MP(1.0))*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/((dn + alpha + beta - MHD_MP(4.0))*dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(1) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(MHD_MP(2.0)*dn + alpha + beta))/(MHD_MP(2.0)*dn);
      cs(2) = ((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0))*(alpha*alpha - beta*beta))/(MHD_MP(2.0)*dn*(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(2.0)));
      cs(3) = MHD_MP(1.0)/(dn + alpha + beta - MHD_MP(3.0));

      cs(0) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))*(dn + alpha + beta - MHD_MP(3.0))*(dn + MHD_MP(1.0))/((MHD_MP(2.0)*dn + alpha + beta - MHD_MP(3.0))*(dn + alpha - MHD_MP(1.0))*(dn + beta - MHD_MP(1.0))));
      cs(1) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(2) *= precision::sqrt((MHD_MP(2.0)*dn + alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*dn + alpha + beta - MHD_MP(1.0)));
      cs(3) *= precision::sqrt((dn + alpha + beta - MHD_MP(2.0))*(dn + MHD_MP(2.0))/((dn + alpha)*(dn + beta)));

      return cs;
   }

   internal::Array WorlandPolynomial::unitWD3P1ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(3);

      cs(0) = (MHD_MP(2.0) + alpha + beta);
      cs(1) = (alpha - beta);
      cs(2) = (alpha + beta + MHD_MP(1.0))/(MHD_MP(2.0)*(alpha + beta - MHD_MP(2.0)));

      cs(2) *= precision::sqrt((alpha + beta + MHD_MP(3.0))*(alpha + beta - MHD_MP(2.0))*(MHD_MP(4.0))/((alpha + beta + MHD_MP(1.0))*(alpha + MHD_MP(1.0))*(beta + MHD_MP(1.0))));
      return cs;
   }

   internal::Array WorlandPolynomial::unitWD3P0ab(const internal::MHDFloat alpha, const internal::MHDFloat beta)
   {
      internal::Array cs(1);

      cs(0) = MHD_MP(0.125)*precision::exp(precisiontr1::lgamma(alpha + beta + MHD_MP(1.0)) - precisiontr1::lgamma(alpha + beta - MHD_MP(2.0)));

      cs(0) *= precision::sqrt(alpha + beta + MHD_MP(1.0))*precision::exp(MHD_MP(0.5)*(precisiontr1::lgamma(alpha + beta - MHD_MP(2.0)) + precisiontr1::lgamma(MHD_MP(4.0)) - precisiontr1::lgamma(alpha + MHD_MP(1.0)) - precisiontr1::lgamma(beta + MHD_MP(1.0))));

      return cs;
   }

   WorlandPolynomial::WorlandPolynomial()
   {
   }

   WorlandPolynomial::~WorlandPolynomial()
   {
   }

}
}
