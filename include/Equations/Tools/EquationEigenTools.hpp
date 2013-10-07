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
#include "Equations/Tools/EquationTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

/**
 * @brief Implementation of some tools for schemes with a single eigen direction
 */
namespace EigenTools {

   /**
    * @brief General implementation of linear row
    */
   template <typename TEquation> static DecoupledZSparse makeLinearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs, const bool hasBoundary);

   /**
    * @brief General implementation of time row
    */
   template <typename TEquation> static DecoupledZSparse makeTimeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs, const bool hasBoundary);

//
// Implementation follows
//

   template <typename TEquation> DecoupledZSparse makeLinearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
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
         SparseMatrix blockMatrix = Tools::makeBlockMatrix(eq.couplingInfo(compId).nBlocks(), eq.couplingInfo(compId).fieldIndex(), colIdx);

         Equations::linearBlock(eq, compId, block, *fIt, eigs, hasBoundary);
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

   template <typename TEquation> DecoupledZSparse makeTimeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      // Storage for the matrix row
      int sysN = eq.couplingInfo(compId).systemN(matIdx);
      DecoupledZSparse  matrixRow(sysN, sysN);
      int blockN = eq.couplingInfo(compId).blockN(matIdx);
      DecoupledZSparse  block(blockN, blockN);
      SparseMatrix  tmp(sysN, sysN);

      // Do only generate equation variable time matrix (for now)
      SpectralFieldId fieldId = std::make_pair(eq.name(), compId);

      // Create time row
      SparseMatrix blockMatrix = Tools::makeBlockMatrix(eq.couplingInfo(compId).nBlocks(), eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex());

      Equations::timeBlock(eq, compId, block, fieldId, eigs, hasBoundary);
      tmp = Eigen::kroneckerProduct(blockMatrix, block.real());
      matrixRow.real() += tmp;
      tmp = Eigen::kroneckerProduct(blockMatrix, block.imag());
      matrixRow.imag() += tmp;

      // Make sure matrices are in compressed format
      matrixRow.real().makeCompressed();
      matrixRow.imag().makeCompressed();

      return matrixRow;
   }

}
}
}

#endif // EQUATIONEIGENTOOLS_HPP
