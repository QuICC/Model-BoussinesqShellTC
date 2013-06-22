/** \file Scalar2DEigenTools.hpp
 *  \brief Implementation of some tools for schemes with two eigen direction
 */

#ifndef SCALAR2DEIGENTOOLS_HPP
#define SCALAR2DEIGENTOOLS_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Project includes
//
#include "Base/Typedefs.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "TypeSelectors/VariableSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of some tools for schemes with two eigen direction
    */
   class Scalar2DEigenTools
   {
      public:
         /**
          * @brief General implementation of linear row for equations with two "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of time row for equations with two "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of boundary row for equations with two "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of the boundary block for equations with two "eigen" dimensions
          */
         static void boundaryBlock2DEigen(const IScalarEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D);
   };

   template <typename TEquation> DecoupledZSparse Scalar2DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = ew.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(sysN, sysN), SparseMatrix(sysN, sysN));
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block = std::make_pair(SparseMatrix(blockN, blockN),SparseMatrix(blockN, blockN));
      SparseMatrix  tmp(sysN, sysN);

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(compId).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
         blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), colIdx) = 1;

         linearBlock(eq, block, *fIt, k2D_, k3D_);
         Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
         matrixRow.first += tmp;
         Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Scalar2DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(sysN, sysN), SparseMatrix(sysN, sysN));
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block = std::make_pair(SparseMatrix(blockN, blockN),SparseMatrix(blockN, blockN));
      SparseMatrix  tmp(sysN, sysN);

      // Create time row
      SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
      blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex()) = 1;
      timeBlock(eq, block, k2D_, 3D_);
      Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
      matrixRow.first += tmp;
      Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Scalar2DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(sysN, sysN), SparseMatrix(sysN, sysN));
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block = std::make_pair(SparseMatrix(blockN, blockN),SparseMatrix(blockN, blockN));
      SparseMatrix  tmp(sysN, sysN);

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(compId).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
         blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), colIdx) = 1;

         boundaryBlock(eq, block, *fIt, k2D_, k3D_);
         Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
         matrixRow.first += tmp;
         Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

}
}

#endif // SCALAR2DEIGENTOOLS_HPP
