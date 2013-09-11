/**
 * @file Equation1DEigenTools.hpp
 * @brief Implementation of some tools for schemes with a single eigen direction 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef EQUATION1DEIGENTOOLS_HPP
#define EQUATION1DEIGENTOOLS_HPP

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
    * @brief Implementation of some tools for schemes with a single eigen direction
    */
   class Equation1DEigenTools
   {
      public:
         /**
          * @brief General implementation of linear row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of time row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of boundary row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx);

         /**
          * @brief General implementation of the boundary block for equations with a single "eigen" dimensions
          */
         static void boundaryBlock1DEigen(const IEquation& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int p1D, const int p3D, const MHDFloat c1D, const MHDFloat c3D);
   };

   template <typename TEquation> DecoupledZSparse Equation1DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

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

         linearBlock(eq, compId, block, *fIt, k_);
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

   template <typename TEquation> DecoupledZSparse Equation1DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(sysN, sysN), SparseMatrix(sysN, sysN));
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block = std::make_pair(SparseMatrix(blockN, blockN),SparseMatrix(blockN, blockN));
      SparseMatrix  tmp(sysN, sysN);

      // Create time row
      SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
      blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex()) = 1;
      timeBlock(eq, compId, block, k_);
      tmp = Eigen::kroneckerProduct(blockMatrix, block.first);
      matrixRow.first += tmp;
      tmp = Eigen::kroneckerProduct(blockMatrix, block.second);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Equation1DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

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

         boundaryBlock(eq, compId, block, *fIt, k_);
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

#endif // EQUATION1DEIGENTOOLS_HPP
