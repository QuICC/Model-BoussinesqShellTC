/**
 * @file EquationEigen2DTools.hpp
 * @brief Implementation of some tools for schemes with two eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGEN2DTOOLS_HPP
#define EQUATIONEIGEN2DTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include<vector>

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/KroneckerTypedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"
#include "TypeSelectors/BoundaryMethodSelector.hpp"
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with two eigen directions
 */
namespace Eigen2D {

   typedef KronEigen2DRProduct KRProduct;
   typedef KronEigen2DZProduct KZProduct;
   typedef KronEigen2DRSum KRSum;
   typedef KronEigen2DZSum KZSum;

   /**
    * @brief Set eigen values
    */
   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx);

   /**
    * @brief Create setup for matrices with minial coupling
    *
    * @param spRes   Shared resolution
    * @param nMat    Number of matrices
    * @param blocks  Size of the blocks
    * @param cols    Number of right-hand sides
    */
   void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols);

   void computeKProduct(SparseMatrix& mat, const KRProduct& block);

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block);

   void computeKSum(SparseMatrix& mat, const KRSum& blocks);

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks);

//
// Implementation follows
//

   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx)
   {
      std::vector<MHDFloat> eigs;

      // Get mode indexes
      ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // k2D_
      eigs.push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0)));
      // k3D_
      eigs.push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1)));

      return eigs;
   }

}
}
}

#endif // EQUATIONEIGEN2DTOOLS_HPP
