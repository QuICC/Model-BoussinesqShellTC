/**
 * @file Equation2DEigenTools.hpp
 * @brief Implementation of some tools for schemes with two eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATION2DEIGENTOOLS_HPP
#define EQUATION2DEIGENTOOLS_HPP

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
#include "Equations/IEquation.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of some tools for schemes with two eigen direction
    */
   class Equation2DEigenTools
   {
      public:
         /**
          * @brief Number of eigen dimensions
          */
         static const int EIGEN_DIMS = 2;

         /**
          * @brief Create setup for matrices with minial coupling
          *
          * @param spRes   Shared resolution
          * @param nMat    Number of matrices
          * @param blocks  Size of the blocks
          * @param cols    Number of right-hand sides
          */
         static void makeMinimalCoupling(const SharedResolution spRes, int& nMat, ArrayI& blocks, ArrayI& cols);

         /**
          * @brief General implementation of linear row
          */
         template <typename TEquation> static DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of time row
          */
         template <typename TEquation> static DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of boundary row
          */
         template <typename TEquation> static DecoupledZSparse boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of the boundary block
          */
         static void boundaryBlock(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat c1D);
   };

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

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

         Equations::linearBlock(eq, compId, block, *fIt, k2D_, k3D_);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
         matrixRow.first += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
         matrixRow.second += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(sysN, sysN), SparseMatrix(sysN, sysN));
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block = std::make_pair(SparseMatrix(blockN, blockN),SparseMatrix(blockN, blockN));
      SparseMatrix  tmp(sysN, sysN);

      // Create time row
      SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
      blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex()) = 1;
      Equations::timeBlock(eq, compId, block, k2D_, k3D_);
      tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
      matrixRow.first += tmp;
      tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Equation2DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get mode indexes
      ArrayI mode = eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(matIdx);

      // Get wave number rescaled to box size (if required)
      MHDFloat k2D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(mode(0));
      MHDFloat k3D_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM3D)*static_cast<MHDFloat>(mode(1));

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

         Equations::boundaryBlock(eq, compId, block, *fIt, k2D_, k3D_);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
         matrixRow.first += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
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

#endif // EQUATION2DEIGENTOOLS_HPP
