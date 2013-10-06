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

   template <typename TEquation> void constrainBlock(TEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, KZSum& blocks, const std::vector<MHDFloat>& bcIdx);

   /**
    * @brief General implementation of linear row
    */
   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief General implementation of time row
    */
   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

   /**
    * @brief General implementation of boundary row
    */
   template <typename TEquation> void boundaryRow(TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);


//   /**
//    * @brief General implementation of the boundary block
//    */
//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D);
//


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

   template <typename TEquation> DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeLinearRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      return EigenTools::makeTimeRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> void boundaryRow(TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      std::vector<MHDFloat> eigs = getEigs(eq, matIdx);

      EigenTools::makeBoundaryRow(eq, compId, matIdx, eigs);
   }

   template <typename TEquation> void storeBoundaryCondition(IEquation& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& coeffs, const std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(coeffs.size() == bcIdx.size());
      assert(coeffs.size() == 1);

      SpectralFieldId eqId = std::make_pair(eq.name(), compId);

      int nEq1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();

      int nI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      Boundary::BCVector bcs1D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second;

      eq.rBcCoord(compId).add1D(bcIdx.at(0), Boundary::MethodSelector<Dimensions::Simulation::SIM1D>::Type(coeffs.at(0), nI, bcs1D, nEq1D));
   }

}
}
}

#endif // EQUATIONEIGEN2DTOOLS_HPP
