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

#if defined GEOMHDISCC_SHNORM_SCHMIDT

   void AssociatedLegendrePolynomial::Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      ipoly.resize(gN, nPoly);
      internal::Array ipmm(gN);
      AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);
      ipoly.col(0) = ipmm;

      if(nPoly > 1)
      {
         internal::Array ipmm1(gN);
         AssociatedLegendrePolynomial::Pmm1(ipmm1, m, ipmm, igrid);
         ipoly.col(1) = ipmm1;
      }

      if(nPoly > 2)
      {
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            ipoly.col(i) = (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipoly.col(i-1).array();
            ipoly.col(i) -= ipoly.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1));
            ipoly.col(i) /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      poly = Precision::cast(ipoly);
   }

   void AssociatedLegendrePolynomial::dPlm(Matrix& diff, internal::Matrix& idiff, const int m, const internal::Matrix& ipoly, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial derivative P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      idiff.resize(gN, nPoly);
      internal::Array idpmm(gN);
      AssociatedLegendrePolynomial::dPmm(idpmm, m, igrid);
      idiff.col(0) = idpmm;

      if(nPoly > 1)
      {
         internal::Array ipmm(gN);
         AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);

         internal::Array idpmm1(gN);
         AssociatedLegendrePolynomial::dPmm1(idpmm1, m, ipmm, idpmm, igrid);
         idiff.col(1) = idpmm1;
      }

      if(nPoly > 2)
      {
         internal::Array isin = igrid.array().acos();
         isin = isin.array().sin();
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            idiff.col(i) = (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*idiff.col(i-1).array();
            idiff.col(i).array() -= (internal::MHDFloat(2*l) - MHD_MP(1.0))*isin.array()*ipoly.col(i-1).array();
            idiff.col(i) -= idiff.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1));
            idiff.col(i) /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      diff = Precision::cast(idiff);
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      ipoly.resize(gN, nPoly);
      internal::Array isin_1pmm(gN);
      AssociatedLegendrePolynomial::sin_1Pmm(isin_1pmm, m, igrid);
      ipoly.col(0) = isin_1pmm;

      if(nPoly > 1)
      {
         internal::Array isin_1pmm1(gN);
         AssociatedLegendrePolynomial::sin_1Pmm1(isin_1pmm1, m, isin_1pmm, igrid);
         ipoly.col(1) = isin_1pmm1;
      }

      if(nPoly > 2)
      {
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            ipoly.col(i) = (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipoly.col(i-1).array();
            ipoly.col(i) -= ipoly.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1));
            ipoly.col(i) /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      poly = Precision::cast(ipoly);
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
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op.setConstant(MHD_MP(0.0));
      } else
      {
         internal::MHDFloat di;
         internal::MHDFloat factor(m);

         for(int i = 1; i <= m; i++)
         {
            di = internal::MHDFloat(2*i);
            factor *= -precision::sqrt((di - MHD_MP(1.0))/di);
         }
         op = igrid.array().acos();
         op = op.array().sin();
         op = op.array().pow(m-1);
         op.array() *= igrid.array();
         op *= factor;
      }
   }

   void AssociatedLegendrePolynomial::dPmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& idpmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid.array().acos();
         op = -op.array().sin();
         op.array() *= ipmm.array();

         op.array() += igrid.array()*idpmm.array();

         op *= precision::sqrt(internal::MHDFloat(2*m + 1));
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op = igrid.array().acos();
         op = op.array().sin().pow(-1);
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
         op = op.array().pow(m-1);
         op *= factor;
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm1(internal::Array& op, const int m, const internal::Array& isin_1pmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op *= precision::sqrt(internal::MHDFloat(2*m + 1));
         op.array() *= isin_1pmm.array(); 
      }
   }

#elif defined GEOMHDISCC_SHNORM_UNITY

   void AssociatedLegendrePolynomial::Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      ipoly.resize(gN, nPoly);
      internal::Array ipmm(gN);
      AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);
      ipoly.col(0) = ipmm;

      if(nPoly > 1)
      {
         internal::Array ipmm1(gN);
         AssociatedLegendrePolynomial::Pmm1(ipmm1, m, ipmm, igrid);
         ipoly.col(1) = ipmm1;
      }

      if(nPoly > 2)
      {
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            ipoly.col(i) = precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipoly.col(i-1).array();
            ipoly.col(i) -= ipoly.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0));
            ipoly.col(i) *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      poly = Precision::cast(ipoly);
   }

   void AssociatedLegendrePolynomial::dPlm(Matrix& diff, internal::Matrix& idiff, const int m, const internal::Matrix& ipoly, const internal::Array& igrid)
   {
      int gN = diff.rows();
      int nPoly = diff.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial derivative P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      idiff.resize(gN, nPoly);
      internal::Array idpmm(gN);
      AssociatedLegendrePolynomial::dPmm(idpmm, m, igrid);
      idiff.col(0) = idpmm;

      if(nPoly > 1)
      {
         internal::Array ipmm(gN);
         AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);

         internal::Array idpmm1(gN);
         AssociatedLegendrePolynomial::dPmm1(idpmm1, m, ipmm, idpmm, igrid);
         idiff.col(1) = idpmm1;
      }

      if(nPoly > 2)
      {
         internal::Array isin = igrid.array().acos();
         isin = isin.array().sin();
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            idiff.col(i) = precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*idiff.col(i-1).array();
            idiff.col(i).array() -= precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*isin.array()*ipoly.col(i-1).array();
            idiff.col(i) -= idiff.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0));

            idiff.col(i).array() *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      diff = Precision::cast(idiff);
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Matrix& poly, internal::Matrix& ipoly, const int m, const internal::Array& igrid)
   {
      int gN = poly.rows();
      int nPoly = poly.cols();

      if (m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      }

      if (nPoly < 1)
      {
         throw Exception("Operator matrix should have at least 1 column");
      }

      ipoly.resize(gN, nPoly);
      internal::Array isin_1pmm(gN);
      AssociatedLegendrePolynomial::sin_1Pmm(isin_1pmm, m, igrid);
      ipoly.col(0) = isin_1pmm;

      if(nPoly > 1)
      {
         internal::Array isin_1pmm1(gN);
         AssociatedLegendrePolynomial::sin_1Pmm1(isin_1pmm1, m, isin_1pmm, igrid);
         ipoly.col(1) = isin_1pmm1;
      }

      if(nPoly > 2)
      {
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            ipoly.col(i) = precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipoly.col(i-1).array();
            ipoly.col(i) -= ipoly.col(i-2)*precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0));
            ipoly.col(i) *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
         }
      }

      poly = Precision::cast(ipoly);
   }

   void AssociatedLegendrePolynomial::Pmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op.setConstant(precision::sqrt(MHD_MP(1.0)/(MHD_MP(4.0)*Precision::PI)));
      } else
      {
         internal::MHDFloat di;
         internal::MHDFloat factor = precision::sqrt((internal::MHDFloat(2*m) + MHD_MP(1.0))/(MHD_MP(4.0)*Precision::PI));

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
         op *= precision::sqrt(internal::MHDFloat(2*m + 3));
         op.array() *= ipmm.array(); 
      }
   }

   void AssociatedLegendrePolynomial::dPmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op.setConstant(MHD_MP(0.0));
      } else
      {
         internal::MHDFloat di;
         internal::MHDFloat factor = internal::MHDFloat(m)*precision::sqrt((internal::MHDFloat(2*m) + MHD_MP(1.0))/(MHD_MP(4.0)*Precision::PI));

         for(int i = 1; i <= m; i++)
         {
            di = internal::MHDFloat(2*i);
            factor *= -precision::sqrt((di - MHD_MP(1.0))/di);
         }
         op = igrid.array().acos();
         op = op.array().sin();
         op = op.array().pow(m-1);
         op.array() *= igrid.array();
         op *= factor;
      }
   }

   void AssociatedLegendrePolynomial::dPmm1(internal::Array& op, const int m, const internal::Array& ipmm, const internal::Array& idpmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid.array().acos();
         op = -op.array().sin();
         op.array() *= ipmm.array();

         op.array() += igrid.array()*idpmm.array();

         op *= precision::sqrt(internal::MHDFloat(2*m + 3));
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm(internal::Array& op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op = igrid.array().acos();
         op = op.array().sin().pow(-1);
         op *= precision::sqrt(MHD_MP(1.0)/(MHD_MP(4.0)*Precision::PI));
      } else
      {
         internal::MHDFloat di;
         internal::MHDFloat factor = precision::sqrt((internal::MHDFloat(2*m) + MHD_MP(1.0))/(MHD_MP(4.0)*Precision::PI));

         for(int i = 1; i <= m; i++)
         {
            di = internal::MHDFloat(2*i);
            factor *= -precision::sqrt((di - MHD_MP(1.0))/di);
         }
         op = igrid.array().acos();
         op = op.array().sin();
         op = op.array().pow(m-1);
         op *= factor;
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm1(internal::Array& op, const int m, const internal::Array& isin_1pmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op *= precision::sqrt(internal::MHDFloat(2*m + 3));
         op.array() *= isin_1pmm.array(); 
      }
   }
#endif //defined GEOMHDISCC_NORMALIZED_SH

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial()
   {
   }

   AssociatedLegendrePolynomial::~AssociatedLegendrePolynomial()
   {
   }

}
}
