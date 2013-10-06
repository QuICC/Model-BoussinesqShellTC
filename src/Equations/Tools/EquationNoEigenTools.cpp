/** 
 * @file EquationNoEigenTools.cpp
 * @brief Source of the tools for schemes with no eigen direction
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
#include "Equations/Tools/EquationNoEigenTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace NoEigen {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);
      // Get 2D dimension (medium)
      int nK = spRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::TRANSFORM);
      // Get 3D dimension (slow)
      int nJ = spRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::TRANSFORM);

      // Get 3D dimension (medium)
      nMat = 1;

      blocks.resize(nMat);
      blocks.setConstant(nI*nJ*nK);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      SparseMatrix tmp = Eigen::kroneckerProduct(std::tr1::get<2>(block), std::tr1::get<0>(block));
      mat = Eigen::kroneckerProduct(std::tr1::get<1>(block), tmp);
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      assert(std::tr1::get<0>(block).imag().nonZeros() == 0);
      assert(std::tr1::get<1>(block).imag().nonZeros() == 0);
      assert(std::tr1::get<2>(block).imag().nonZeros() == 0);

      SparseMatrix tmp = Eigen::kroneckerProduct(std::tr1::get<2>(block).real(), std::tr1::get<0>(block).real());
      mat.real() = Eigen::kroneckerProduct(std::tr1::get<1>(block).real(), tmp);
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

//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D)
//   {
//      throw Exception("Not yet implemented!");
//   }
}
}
}
