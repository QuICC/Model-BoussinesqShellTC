/**
 * @file EquationNoEigenTools.hpp
 * @brief Implementation of some tools for schemes with no eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONNOEIGENTOOLS_HPP
#define EQUATIONNOEIGENTOOLS_HPP

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
 * @brief Tools for equations with no eigen direction
 */
namespace NoEigen {

   typedef KronNoEigenRProduct KRProduct;
   typedef KronNoEigenZProduct KZProduct;
   typedef KronNoEigenRSum KRSum;
   typedef KronNoEigenZSum KZSum;

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
      throw Exception("Not yet implemented!");
      std::vector<MHDFloat> eigs;

      return eigs;
   }

}
}

#endif // EQUATIONNOEIGENTOOLS_HPP
