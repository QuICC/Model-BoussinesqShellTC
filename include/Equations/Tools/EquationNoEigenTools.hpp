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
#include "Equations/Tools/EquationEigenTools.hpp"

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

//
//   /**
//    * @brief General implementation of the boundary block
//    */
//   void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D);


//
// Implementation follows
//

   template <typename TEquation> std::vector<MHDFloat> getEigs(const TEquation& eq, const int matIdx)
   {
      throw Exception("Not yet implemented!");
      std::vector<MHDFloat> eigs;

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
      assert(coeffs.size() == 3);

      SpectralFieldId eqId = std::make_pair(eq.name(), compId);

      int nEq1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();
      int nEq2D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM2D)->second.size();
      int nEq3D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM3D)->second.size();

      int nI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::SPECTRAL);
      int nJ = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      Boundary::BCVector bcs1D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second;
      Boundary::BCVector bcs2D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM2D)->second;
      Boundary::BCVector bcs3D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM3D)->second;

      eq.rBcCoord(compId).add1D(bcIdx.at(0), Boundary::MethodSelector<Dimensions::Simulation::SIM1D>::Type(coeffs.at(0), nI, bcs1D, nEq1D));
      eq.rBcCoord(compId).add2D(bcIdx.at(1), Boundary::MethodSelector<Dimensions::Simulation::SIM2D>::Type(coeffs.at(1), nI, bcs2D, nEq2D));
      eq.rBcCoord(compId).add3D(bcIdx.at(2), Boundary::MethodSelector<Dimensions::Simulation::SIM3D>::Type(coeffs.at(2), nI, bcs3D, nEq3D));
   }

}
}

#endif // EQUATIONNOEIGENTOOLS_HPP
