/**
 * @file EquationEigenTools.hpp
 * @brief Implementation of generic tools to deal with eigen dimensions in equations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGENTOOLS_HPP
#define EQUATIONEIGENTOOLS_HPP

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
#include "Equations/EquationTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   /**
    * @brief Implementation of some tools for schemes with a single eigen direction
    */
   class EquationEigenTools
   {
      public:
         /**
          * @brief General implementation of linear row
          */
         template <typename TEquation> static DecoupledZSparse makeLinearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs);

         /**
          * @brief General implementation of time row
          */
         template <typename TEquation> static DecoupledZSparse makeTimeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs);

         /**
          * @brief General implementation of boundary row
          */
         template <typename TEquation> static DecoupledZSparse makeBoundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs);
   };

   template <typename TEquation> DecoupledZSparse EquationEigenTools::makeLinearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs)
   {
      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow(sysN, sysN);
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block(blockN, blockN);
      SparseMatrix  tmp(sysN, sysN);

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(compId).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix blockMatrix = internal::makeBlockMatrix(eq.couplingInfo(compId).nBlocks(), eq.couplingInfo(compId).fieldIndex(), colIdx);

         Equations::linearBlock(eq, compId, block, *fIt, eigs);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.real());
         matrixRow.real() += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.imag());
         matrixRow.imag() += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.real().makeCompressed();
      matrixRow.imag().makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse EquationEigenTools::makeTimeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs)
   {
      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow(sysN, sysN);
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block(blockN, blockN);
      SparseMatrix  tmp(sysN, sysN);

      // Create time row
      SparseMatrix blockMatrix = internal::makeBlockMatrix(eq.couplingInfo(compId).nBlocks(), eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex());

      Equations::timeBlock(eq, compId, block, eigs);
      tmp = Eigen::kroneckerProduct(blockMatrix, block.real());
      matrixRow.real() += tmp;
      tmp = Eigen::kroneckerProduct(blockMatrix, block.imag());
      matrixRow.imag() += tmp;

      // Make sure matrices are in compressed format
      matrixRow.real().makeCompressed();
      matrixRow.imag().makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse EquationEigenTools::makeBoundaryRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs)
   {
      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow(sysN, sysN);
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block(blockN, blockN);
      SparseMatrix  tmp(sysN, sysN);

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(compId).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix blockMatrix = internal::makeBlockMatrix(eq.couplingInfo(compId).nBlocks(), eq.couplingInfo(compId).fieldIndex(), colIdx);

         Equations::boundaryBlock(eq, compId, block, *fIt, eigs);
         tmp = Eigen::kroneckerProduct(blockMatrix, block.real());
         matrixRow.real() += tmp;
         tmp = Eigen::kroneckerProduct(blockMatrix, block.imag());
         matrixRow.imag() += tmp;

         colIdx++;
      }

      // Make sure matrices are in compressed format
      matrixRow.real().makeCompressed();
      matrixRow.imag().makeCompressed();

      return matrixRow;
   }

}
}

#endif // EQUATIONEIGENTOOLS_HPP
