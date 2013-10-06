/** 
 * @file EquationEigen3DTools.cpp
 * @brief Source of the tools for schemes with three eigen direction
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
#include "Equations/Tools/EquationEigen3DTools.hpp"

// Project includes
//
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

namespace Eigen3D {

   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols)
   {
      // Get 1D dimension (fast)
      int nI = spRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::TRANSFORM);

      for(int i = 0; i < spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         nMat += spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i)*spRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATB1D>(i,j);
      }

      blocks.resize(nMat);
      blocks.setConstant(1);
      cols.resize(nMat);
      cols.setConstant(1);
   }

   void computeKProduct(SparseMatrix& mat, const KRProduct& block)
   {
      mat.resize(1,1);
      mat(0,0) = block;
   }

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block)
   {
      mat.resize(1,1);
      mat(0,0).real() = block.real();
      mat(0,0).imag() = block.imag();
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
