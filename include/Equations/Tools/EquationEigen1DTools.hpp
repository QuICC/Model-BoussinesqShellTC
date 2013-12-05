/**
 * @file EquationEigen1DTools.hpp
 * @brief Implementation of some tools for schemes with a single eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGEN1DTOOLS_HPP
#define EQUATIONEIGEN1DTOOLS_HPP

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
#include "Equations/Tools/EquationEigenTools.hpp"
#include "SpectralOperators/UnitOperator.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Tools for equations with a single eigen direction
 */
namespace Eigen1D {

   typedef KronEigen1DRProduct KRProduct;
   typedef KronEigen1DZProduct KZProduct;
   typedef KronEigen1DRSum KRSum;
   typedef KronEigen1DZSum KZSum;

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

   void makeComplex(KZSum &zBlock, const KRSum& rBlock);

   void computeKProduct(SparseMatrix& mat, const KRProduct& block);

   void computeKProduct(DecoupledZSparse& mat, const KZProduct& block);

   void computeKSum(SparseMatrix& mat, const KRSum& blocks);

   void computeKSum(DecoupledZSparse& mat, const KZSum& blocks);

   template <typename TEquation> void constrainBlock(TEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, KZSum& blocks, const std::vector<MHDFloat>& bcIdx, const int hasBoundary);

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

      // Get wave number rescale to box size
      eigs.push_back(eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx)));

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

   template <typename TEquation> void constrainBlock(TEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, KZSum& blocks, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      // Get boundary information
      std::vector<MHDFloat> coeffs;
      std::vector<Boundary::BCIndex> bcIdx;
      Equations::boundaryBlock(eq, compId, fieldId, eigs, coeffs, bcIdx);

      assert(coeffs.size() == bcIdx.size());
      assert(coeffs.size() == 2);

      SpectralFieldId eqId = std::make_pair(eq.name(), compId);
      int nEq1D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM1D)->second.size();
      int nEq3D = eq.bcIds().bcs(eqId,eqId).find(Dimensions::Simulation::SIM3D)->second.size();

      int nI = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nK = eq.spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      Boundary::BCVector bcs1D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM1D)->second;
      Boundary::BCVector bcs3D = eq.bcIds().bcs(eqId,fieldId).find(Dimensions::Simulation::SIM3D)->second;

      Boundary::MethodSelector<Dimensions::Simulation::SIM1D>::Type bcOp1D(coeffs.at(0), nI, bcs1D, nEq1D);
      Boundary::MethodSelector<Dimensions::Simulation::SIM3D>::Type bcOp3D(coeffs.at(1), nK, bcs3D, nEq3D);

      for(KZSum::iterator it = blocks.begin(); it != blocks.end(); ++it)
      {
         bcOp1D.constrainKronBlock(std::tr1::get<0>(*it));
         bcOp3D.constrainKronBlock(std::tr1::get<1>(*it));
      }

      if(hasBoundary)
      {
         DecoupledZSparse tmp;
         Spectral::UnitOperator idOp(bcOp3D.nN());
         if(bcOp1D.constrainKronProduct(tmp))
         {
            DecoupledZSparse idMat(bcOp3D.nN(),bcOp3D.nN());
            idMat.real() = idOp.id();
            blocks.push_back(std::tr1::make_tuple(tmp, idMat));
         }
         if(bcOp3D.constrainKronProduct(tmp))
         {
            idOp.reset(bcOp1D.nN());
            DecoupledZSparse idMat(bcOp1D.nN(),bcOp1D.nN());
            idMat.real() = idOp.shiftId(bcOp3D.nBc());
            blocks.push_back(std::tr1::make_tuple(idMat, tmp));
         }
      }

      computeKSum(mat, blocks);
   }

}
}
}

#endif // EQUATIONEIGEN1DTOOLS_HPP
