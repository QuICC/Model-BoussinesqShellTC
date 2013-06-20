/** \file Scalar1DEigenTools.hpp
 *  \brief Implementation of some tools for schemes with a single eigen direction
 */

#ifndef SCALAR1DEIGENTOOLS_HPP
#define SCALAR1DEIGENTOOLS_HPP

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
    * @brief Implementation of some tools for schemes with a single eigen direction
    */
   class Scalar1DEigenTools
   {
      public:
         /**
          * @brief General implementation of linear row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse linearRow(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

         /**
          * @brief General implementation of time row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse timeRow(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

         /**
          * @brief General implementation of boundary row for equations with a single "eigen" dimension
          */
         template <typename TEquation> static DecoupledZSparse boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx);

         /**
          * @brief General implementation of the boundary block for equations with a single "eigen" dimensions
          */
         static void boundaryBlock(const Scalar1DEigenTools& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const int pX, const int pZ, const MHDFloat cX, const MHDFloat cZ);
   };

   template <typename TEquation> DecoupledZSparse Scalar1DEigenTools::linearRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)), SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)),SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
         blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), colIdx) = 1;

         linearBlock(eq, block, *fIt, k_);
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

   template <typename TEquation> DecoupledZSparse Scalar1DEigenTools::timeRow(const TEquation& eq, FieldComponents::Spectral::Id compId, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)), SparseMatrix(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)),SparseMatrix(eq.couplingInfo(compId).blockN(matIdx), eq.couplingInfo(compId).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(compId).systemN(matIdx), eq.couplingInfo(compId).systemN(matIdx));

      // Create time row
      SparseMatrix   blockMatrix(eq.couplingInfo(compId).nBlocks(),eq.couplingInfo(compId).nBlocks());
      blockMatrix.insert(eq.couplingInfo(compId).fieldIndex(), eq.couplingInfo(compId).fieldIndex()) = 1;
      timeBlock(eq, block, k_);
      Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
      matrixRow.first += tmp;
      Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   template <typename TEquation> DecoupledZSparse Scalar1DEigenTools::boundaryRow(const TEquation& eq, FieldComponents::Spectral::Id comp, const int matIdx)
   {
      // Get wave number rescale to box size
      MHDFloat k_ = eq.unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(eq.unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->template idx<Dimensions::Data::DAT3D>(matIdx));

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx)), SparseMatrix(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(eq.couplingInfo(comp).blockN(matIdx), eq.couplingInfo(comp).blockN(matIdx)),SparseMatrix(eq.couplingInfo(comp).blockN(matIdx), eq.couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(eq.couplingInfo(comp).systemN(matIdx), eq.couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::FieldId_iterator fIt;
      CouplingInformation::FieldId_range fRange = eq.couplingInfo(FieldComponents::Spectral::SCALAR).implicitRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(eq.couplingInfo(comp).nBlocks(),eq.couplingInfo(comp).nBlocks());
         blockMatrix.insert(eq.couplingInfo(comp).fieldIndex(), colIdx) = 1;

         boundaryBlock(eq, block, *fIt, k_);
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

#endif // SCALAR1DEIGENTOOLS_HPP
