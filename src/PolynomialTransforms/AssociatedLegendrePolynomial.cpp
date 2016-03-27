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
         AssociatedLegendrePolynomial::Pm1m(ipoly.col(1), m, ipoly.col(0), igrid);
      }

      if(nPoly > 2)
      {
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            AssociatedLegendrePolynomial::Plm(ipoly.col(i), m, l, ipoly.col(i-1), ipoly.col(i-2), igrid);
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
         AssociatedLegendrePolynomial::dPm1m(idiff.col(1), m, ipoly.col(0), idiff.col(0), igrid);
      }

      if(nPoly > 2)
      {
         internal::Array isin = igrid.array().acos();
         isin = isin.array().sin();
         for(int i = 2; i < nPoly; ++i)
         {
            int l = m + i;
            AssociatedLegendrePolynomial::dPlm(idiff.col(i), m, l, idiff.col(i-1), idiff.col(i-2), ipoly.col(i-1), igrid);
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

      // Polynomials is set to zero for m=0 as it only appears combined with \partial_\phi
      if(m == 0)
      {
         ipoly.resize(gN, nPoly);
         ipoly.setZero();

         poly.resize(gN, nPoly);
         poly.setZero();

      } else
      {
         // Storage for P_{l+1}^{m+1} and P_{l+1}^{m-1}
         internal::Matrix ipl1m1(gN, 2);
         internal::Matrix ipl1m_1(gN, 2);

         // Initialize P_{l+1}^{m+1}
         AssociatedLegendrePolynomial::Pmm(ipl1m1.col(0), m+1, igrid);

         // Initialize P_{l+1}^{m-1}
         AssociatedLegendrePolynomial::Pmm(ipl1m_1.col(0), m-1, igrid);
         AssociatedLegendrePolynomial::Pm1m(ipl1m_1.col(1), m-1, ipl1m_1.col(0), igrid);
         AssociatedLegendrePolynomial::Plm(ipl1m_1.col(0), m-1, m+1, ipl1m_1.col(1), ipl1m_1.col(0), igrid);
         ipl1m_1.col(0).swap(ipl1m_1.col(1));

         // Initialize \frac{1}{sin(\theta)} P_{l}^{m}
         ipoly.resize(gN, nPoly);
         AssociatedLegendrePolynomial::sin_1Plm(ipoly.col(0), m, m, ipl1m1.col(0), ipl1m_1.col(1));

         if(nPoly > 1)
         {
            // Increment P_{l+1}^{m+1}
            AssociatedLegendrePolynomial::Pm1m(ipl1m1.col(1), m+1, ipl1m1.col(0), igrid);

            // Increment P_{l+1}^{m-1}
            AssociatedLegendrePolynomial::Plm(ipl1m_1.col(0), m-1, m+2, ipl1m_1.col(1), ipl1m_1.col(0), igrid);
            ipl1m_1.col(0).swap(ipl1m_1.col(1));

            // Increment \frac{1}{sin(\theta)} P_{l}^{m}
            AssociatedLegendrePolynomial::sin_1Plm(ipoly.col(1), m, m+1, ipl1m1.col(1), ipl1m_1.col(1));
         }

         if(nPoly > 2)
         {
            for(int i = 2; i < nPoly; ++i)
            {
               int l = m + i;

               // Increment P_{l+1}^{m+1}
               AssociatedLegendrePolynomial::Plm(ipl1m1.col(0), m+1, l+1, ipl1m1.col(1), ipl1m1.col(0), igrid);
               ipl1m1.col(0).swap(ipl1m1.col(1));

               // Increment P_{l+1}^{m-1}
               AssociatedLegendrePolynomial::Plm(ipl1m_1.col(0), m-1, l+1, ipl1m_1.col(1), ipl1m_1.col(0), igrid);
               ipl1m_1.col(0).swap(ipl1m_1.col(1));

               // Increment \frac{1}{sin(\theta)} P_{l}^{m}
               AssociatedLegendrePolynomial::sin_1Plm(ipoly.col(i), m, l, ipl1m1.col(1), ipl1m_1.col(1));
            }
         }

         poly = Precision::cast(ipoly);
      }
   }

#if defined GEOMHDISCC_SHNORM_SCHMIDT

   void AssociatedLegendrePolynomial::Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);

      iplm.array() = -precision::sqrt((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(internal::MHDFloat(2*l) - MHD_MP(1.0))*ipl_2m.array();
      iplm.array() += igrid.array()*ipl_1m.array();
      iplm.array() *= (internal::MHDFloat(2*l) - MHD_MP(1.0))/precision::sqrt(dl*dl - dm*dm);
   }

   void AssociatedLegendrePolynomial::dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& idpl_1m, const Eigen::Ref<const internal::Matrix>& idpl_2m, const Eigen::Ref<const internal::Matrix>& ipl_1m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::Array isin = igrid.array().acos();
      isin = isin.array().sin();

      idplm.array() = -precision::sqrt(internal::MHDFloat(l + m - 1)*internal::MHDFloat(l - m - 1))/(internal::MHDFloat(2*l) - MHD_MP(1.0))*idpl_2m.array();
      idplm.array() += igrid.array()*idpl_1m.array();
      idplm.array() -= isin.array()*ipl_1m.array();
      idplm.array() *= (internal::MHDFloat(2*l) - MHD_MP(1.0))/precision::sqrt(dl*dl - dm*dm);
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl1m1, const Eigen::Ref<const internal::Matrix>& ipl1m_1)
   {
      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);

      iplm.array() = ipl1m1.array();
      iplm.array() += precision::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))))*ipl1m_1.array();
      iplm.array() *= -precision::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/MHD_MP(2.0*dm);
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

   void AssociatedLegendrePolynomial::Pm1m(Eigen::Ref<internal::Matrix> op, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid)
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

   void AssociatedLegendrePolynomial::dPm1m(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& ipmm, const Eigen::Ref<const internal::Array>& idpmm, const internal::Array& igrid)
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

#elif defined GEOMHDISCC_SHNORM_UNITY

   void AssociatedLegendrePolynomial::Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl_1m, const Eigen::Ref<const internal::Matrix>& ipl_2m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);

      iplm.array() = -precision::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)))*ipl_2m.array();
      iplm.array() += igrid.array()*ipl_1m.array();
      iplm.array() *= precision::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));
   }

   void AssociatedLegendrePolynomial::dPlm(Eigen::Ref<internal::Matrix> idplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& idpl_1m, const Eigen::Ref<const internal::Matrix>& idpl_2m, const Eigen::Ref<const internal::Matrix>& ipl_1m, const internal::Array& igrid)
   {
      // Safety assert
      assert(l-m > 1);

      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);

      internal::Array isin = igrid.array().acos();
      isin = isin.array().sin();

      idplm.array() = -precision::sqrt(((dl - MHD_MP(1.0))*(dl - MHD_MP(1.0)) - dm*dm)/(MHD_MP(4.0)*dl*(dl - MHD_MP(2.0)) + MHD_MP(3.0)))*idpl_2m.array();
      idplm.array() += igrid.array()*idpl_1m.array();
      idplm.array() -= isin.array()*ipl_1m.array();
      idplm.array() *= precision::sqrt((MHD_MP(4.0)*dl*dl - MHD_MP(1.0))/(dl*dl - dm*dm));
   }

   void AssociatedLegendrePolynomial::sin_1Plm(Eigen::Ref<internal::Matrix> iplm, const int m, const int l, const Eigen::Ref<const internal::Matrix>& ipl1m1, const Eigen::Ref<const internal::Matrix>& ipl1m_1)
   {
      internal::MHDFloat dl = internal::MHDFloat(l);
      internal::MHDFloat dm = internal::MHDFloat(m);

      iplm.array() = ipl1m1.array();
      iplm.array() += precision::sqrt(((dl - dm + MHD_MP(1.0))*(dl - dm + MHD_MP(2.0)))/((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0))))*ipl1m_1.array();
      iplm.array() *= -precision::sqrt((MHD_MP(2.0)*dl + MHD_MP(1.0))/(MHD_MP(2.0)*dl + MHD_MP(3.0)))*precision::sqrt((dl + dm + MHD_MP(1.0))*(dl + dm + MHD_MP(2.0)))/(2.0*dm);
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

   void AssociatedLegendrePolynomial::Pm1m(Eigen::Ref<internal::Matrix> op, const int m, const Eigen::Ref<const internal::Matrix>& ipmm, const internal::Array& igrid)
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

   void AssociatedLegendrePolynomial::dPm1m(Eigen::Ref<internal::Array> op, const int m, const Eigen::Ref<const internal::Array>& ipmm, const Eigen::Ref<const internal::Array>& idpmm, const internal::Array& igrid)
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
#endif //defined GEOMHDISCC_SHNORM_SCHMIDT

   AssociatedLegendrePolynomial::AssociatedLegendrePolynomial()
   {
   }

   AssociatedLegendrePolynomial::~AssociatedLegendrePolynomial()
   {
   }

}
}
