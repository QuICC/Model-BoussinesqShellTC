/** 
 * @file EquationEigen2DTools.cpp
 * @brief Source of the tools for schemes with two eigen direction
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
#include "Equations/Tools/EquationEigen2DTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen2D {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      blocks.resize(nMat);
      blocks.setConstant(nI);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   int fieldCouplingNMat(const SharedResolution spRes)
   {
      int nMat = 0;

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
      }

      return nMat;
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      assert(block.size() > 0);

      mat = block;
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      assert(block.real().size() > 0);
      assert(block.real().size() == block.imag().size());

      mat.real() = block.real();
      mat.imag() = block.imag();
   }

   void computeKSum(SparseMatrix& mat, const KRSum& blocks)
   {
      computeKProduct(mat, blocks);
   }

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks)
   {
      computeKProduct(mat, blocks);
   }

}
}
}
