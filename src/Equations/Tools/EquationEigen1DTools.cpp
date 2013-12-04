/** 
 * @file EquationEigen1DTools.cpp
 * @brief Source of the tools for schemes with a single eigen direction
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <limits>

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Tools/EquationEigen1DTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen1D {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);
      // Get 2D dimension (slow)
      int nJ = spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::TRANSFORM);

      // Get 3D dimension (medium)
      nMat = spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

      blocks.resize(nMat);
      blocks.setConstant(nI*nJ);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void makeComplex(KZSum& zBlock, const KRSum& rBlock)
   {
      DecoupledZSparse  bL;
      DecoupledZSparse  bR;

      for(KRSum::const_iterator it = rBlock.begin(); it != rBlock.end(); ++it)
      {
         bL.real() = std::tr1::get<0>(*it);
         bR.real() = std::tr1::get<1>(*it);
         zBlock.push_back(std::tr1::make_tuple(bL, bR));
      }
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      assert(std::tr1::get<0>(block).size() > 0);
      assert(std::tr1::get<1>(block).size() > 0);

      mat = Eigen::kroneckerProduct(std::tr1::get<1>(block), std::tr1::get<0>(block));
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      assert(std::tr1::get<0>(block).real().size() > 0);
      assert(std::tr1::get<1>(block).real().size() > 0);
      assert(std::tr1::get<0>(block).real().size() == std::tr1::get<0>(block).imag().size());
      assert(std::tr1::get<1>(block).real().size() == std::tr1::get<1>(block).imag().size());

      mat.real() = Eigen::kroneckerProduct(std::tr1::get<1>(block).real(), std::tr1::get<0>(block).real());
      mat.real() -= Eigen::kroneckerProduct(std::tr1::get<1>(block).imag(), std::tr1::get<0>(block).imag());
      mat.imag() = Eigen::kroneckerProduct(std::tr1::get<1>(block).real(), std::tr1::get<0>(block).imag());
      mat.imag() += Eigen::kroneckerProduct(std::tr1::get<1>(block).imag(), std::tr1::get<0>(block).real());
   }

   void computeKSum(SparseMatrix& mat, const KRSum& blocks)
   {
      if(blocks.size() > 0)
      {
         computeKProduct(mat, blocks.at(0));

         SparseMatrix tmp;
         for(std::vector<KRProduct>::const_iterator it = (blocks.begin()+1); it != blocks.end(); ++it)
         {
            computeKProduct(tmp, *it);

            mat += tmp;
         }
      }
   }

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks)
   {
      if(blocks.size() > 0)
      {
         computeKProduct(mat, blocks.at(0));

         DecoupledZSparse tmp;
         for(std::vector<KZProduct>::const_iterator it = (blocks.begin()+1); it != blocks.end(); ++it)
         {
            computeKProduct(tmp, *it);

            mat.real() += tmp.real();
            mat.imag() += tmp.imag();
         }
      }
   }
}
}
}
