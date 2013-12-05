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
#include "Equations/Tools/EquationEigenTools.hpp"

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

   template <typename TEquation> void constrainBlock(TEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, KZSum& blocks, const std::vector<MHDFloat>& eigs, const int hasBoundary);

   /**
    * @brief General implementation of linear row
    */
   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary);

   /**
    * @brief General implementation of time row
    */
   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary);

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

   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeLinearRow(eq, compId, matIdx, eigs, hasBoundary);
   }

   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeTimeRow(eq, compId, matIdx, eigs, hasBoundary);
   }

   template <typename TEquation> void constrainBlock(TEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, KZSum& blocks, const std::vector<MHDFloat>& eigs, const int hasBoundary)
   {
      // Get boundary information
      std::vector<MHDFloat> coeffs;
      std::vector<Boundary::BCIndex> bcIdx;
      Equations::boundaryBlock(eq, compId, fieldId, eigs, coeffs, bcIdx);

      assert(coeffs.size() == bcIdx.size());
      assert(coeffs.size() == 1);

      SpectralFieldId eqId = std::make_pair(eq.name(), compId);
      int nEq1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();

      int nI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      Boundary::BCVector bcs1D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second;

      Boundary::MethodSelector<Dimensions::Simulation::SIM1D>::Type bcOp1D(coeffs.at(0), nI, bcs1D, nEq1D);

      bcOp1D.constrainKronBlock(blocks);

      if(hasBoundary)
      {
      }

      computeKSum(mat, blocks);
   }

}
}
}

#endif // EQUATIONEIGEN2DTOOLS_HPP
