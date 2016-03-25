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
      AssociatedLegendrePolynomial::Pmm(ipoly.col(0), m, igrid);

      if(nPoly > 1)
      {
         AssociatedLegendrePolynomial::Pmm1(ipoly.col(1), m, ipoly.col(0), igrid);
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
      AssociatedLegendrePolynomial::dPmm(idiff.col(0), m, igrid);

      if(nPoly > 1)
      {
         internal::Array ipmm(gN);
         AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);

         AssociatedLegendrePolynomial::dPmm1(idiff.col(1), m, ipmm, idiff.col(0), igrid);
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
      AssociatedLegendrePolynomial::sin_1Pmm(ipoly.col(0), m, igrid);

      if(nPoly > 1)
      {
         AssociatedLegendrePolynomial::sin_1Pmm1(ipoly.col(1), m, ipoly.col(0), igrid);
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

   void AssociatedLegendrePolynomial::Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      iplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))*ipl_2m.array();
      iplm.array() += (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipl_1m.array();
      iplm.array() /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& idpl_1m, const Eigen::Ref<const internal::Matrix>& idpl_2m, const Eigen::Ref<const internal::Matrix>& ipl_1m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::Array isin = igrid.array().acos();
      isin = isin.array().sin();

      idplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))*idpl_2m.array();
      idplm.array() += (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*idpl_1m.array();
      idplm.array() -= (internal::MHDFloat(2*l) - MHD_MP(1.0))*isin.array()*ipl_1m.array();
      idplm.array() /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      iplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))*ipl_2m.array();
      iplm.array() += (internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipl_1m.array();
      iplm.array() /= precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::Pmm(Eigen::Ref<internal::Matrix> op, const int m, const internal::Array& igrid)
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::Pmm1(Eigen::Ref<internal::Matrix> op, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 1));
         op.array() *= ipmm.array(); 
      }
   }

   void AssociatedLegendrePolynomial::dPmm(Eigen::Ref<internal::Array> op, const int m, const internal::Array& igrid)
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::dPmm1(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& ipmm, const Eigen::Ref<const internal::Array>& idpmm, const internal::Array& igrid)
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

         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 1));
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm(Eigen::Ref<internal::Array> op, const int m, const internal::Array& igrid)
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm1(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& isin_1pmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 1));
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
      AssociatedLegendrePolynomial::Pmm(ipoly.col(0), m, igrid);

      if(nPoly > 1)
      {
         AssociatedLegendrePolynomial::Pmm1(ipoly.col(1), m, ipoly.col(0), igrid);
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
      AssociatedLegendrePolynomial::dPmm(idiff.col(0), m, igrid);

      if(nPoly > 1)
      {
         internal::Array ipmm(gN);
         AssociatedLegendrePolynomial::Pmm(ipmm, m, igrid);

         AssociatedLegendrePolynomial::dPmm1(idiff.col(1), m, ipmm, idiff.col(0), igrid);
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
      AssociatedLegendrePolynomial::sin_1Pmm(ipoly.col(0), m, igrid);

      if(nPoly > 1)
      {
         AssociatedLegendrePolynomial::sin_1Pmm1(ipoly.col(1), m, ipoly.col(0), igrid);
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

   void AssociatedLegendrePolynomial::Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      iplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0))*ipl_2m.array();
      iplm.array() += precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipl_1m.array();
      iplm.array() *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& idpl_1m, const Eigen::Ref<const internal::Matrix>& idpl_2m, const Eigen::Ref<const internal::Matrix>& ipl_1m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::Array isin = igrid.array().acos();
      isin = isin.array().sin();

      idplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0))*idpl_2m.array();
      idplm.array() += precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*idpl_1m.array();
      idplm.array() -= precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*isin.array()*ipl_1m.array();

      idplm.array() *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      iplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(3.0))*ipl_2m.array();
      iplm.array() += precision::sqrt(internal::MHDFloat(2*l) - MHD_MP(1.0))*igrid.array()*ipl_1m.array();
      iplm.array() *= precision::sqrt(internal::MHDFloat(2*l) + MHD_MP(1.0))/precision::sqrt(internal::MHDFloat(l + m)*internal::MHDFloat(l - m));
   }

   void AssociatedLegendrePolynomial::Pmm(Eigen::Ref<internal::Matrix> op, const int m, const internal::Array& igrid)
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::Pmm1(Eigen::Ref<internal::Matrix> op, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 3));
         op.array() *= ipmm.array(); 
      }
   }

   void AssociatedLegendrePolynomial::dPmm(Eigen::Ref<internal::Array> op, const int m, const internal::Array& igrid)
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::dPmm1(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& ipmm, const Eigen::Ref<const internal::Array>& idpmm, const internal::Array& igrid)
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

         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 3));
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm(Eigen::Ref<internal::Array> op, const int m, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else if(m == 0)
      {
         op = igrid.array().acos();
         op = op.array().sin().pow(-1);
         op.array() *= precision::sqrt(MHD_MP(1.0)/(MHD_MP(4.0)*Precision::PI));
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
         op.array() *= factor;
      }
   }

   void AssociatedLegendrePolynomial::sin_1Pmm1(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& isin_1pmm, const internal::Array& igrid)
   {
      if(m < 0)
      {
         throw Exception("Tried to compute associated Legendre polynomial P_l^m with m < 0");
      } else
      {
         op = igrid;
         op.array() *= precision::sqrt(internal::MHDFloat(2*m + 3));
         op.array() *= isin_1pmm.array(); 
      }
   }
#endif //defined GEOMHDISCC_SHNORM_SCHMIDT

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial()
   {
   }

   AssociatedLegendrePolynomial::~AssociatedLegendrePolynomial()
   {
   }

}
}
