/** \file IBeta3DQGScalarEquation.cpp
 *  \brief Source of the implementation of the streamfunction equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/IBeta3DQGScalarEquation.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"

#include <iostream>

namespace GeoMHDiSCC {

namespace Equations {

   IBeta3DQGScalarEquation::IBeta3DQGScalarEquation(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   IBeta3DQGScalarEquation::~IBeta3DQGScalarEquation()
   {
   }

   void IBeta3DQGScalarEquation::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Set coupling information
      this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      Beta3DQGSystem::setCouplingInfo(this->mCouplingInfos.find(FieldComponents::Spectral::SCALAR)->second, this->name(), nX, nZ, nY);
   }

   DecoupledZSparse IBeta3DQGScalarEquation::linearRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get X and Z dimensions
      int nx = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nz = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get wave number rescale to box size
      MHDFloat k = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::field_iterator fIt;
      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).fieldRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
         SparseMatrix   blockMatrix(this->couplingInfo(comp).nBlocks(),this->couplingInfo(comp).nBlocks());
         blockMatrix.insert(this->couplingInfo(comp).fieldIndex(), colIdx) = 1;

         Beta3DQGSystem::linearBlock(block, this->name(), fIt->first, nx, nz, k, Ra, Pr, Gamma, chi);
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

   DecoupledZSparse IBeta3DQGScalarEquation::timeRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get X and Z dimensions
      int nx = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nz = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get wave number rescale to box size
      MHDFloat k = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Create time row
      SparseMatrix   blockMatrix(this->couplingInfo(comp).nBlocks(),this->couplingInfo(comp).nBlocks());
      blockMatrix.insert(this->couplingInfo(comp).fieldIndex(), this->couplingInfo(comp).fieldIndex()) = 1;
      Beta3DQGSystem::timeBlock(block, this->name(), nx, nz, k, Ra, Pr, Gamma, chi);
      Eigen::kroneckerProduct(blockMatrix, block.first, tmp);
      matrixRow.first += tmp;
      Eigen::kroneckerProduct(blockMatrix, block.second, tmp);
      matrixRow.second += tmp;

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   DecoupledZSparse IBeta3DQGScalarEquation::boundaryRow(FieldComponents::Spectral::Id comp, const int matIdx) const
   {
      // Get X and Z dimensions
      int nx = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nz = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Get wave number rescale to box size
      MHDFloat k = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D)*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = this->eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = this->eqParams().nd(NonDimensional::PRANDTL);
      MHDFloat Gamma = this->eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = this->eqParams().nd(NonDimensional::CHI);

      // Storage for the matrix row
      DecoupledZSparse  matrixRow = std::make_pair(SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)), SparseMatrix(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx)));
      DecoupledZSparse  block = std::make_pair(SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)),SparseMatrix(this->couplingInfo(comp).blockN(matIdx), this->couplingInfo(comp).blockN(matIdx)));
      SparseMatrix  tmp(this->couplingInfo(comp).systemN(matIdx), this->couplingInfo(comp).systemN(matIdx));

      // Loop over all coupled fields
      int colIdx = 0;
      CouplingInformation::field_iterator fIt;
      CouplingInformation::field_iterator_range fRange = this->couplingInfo(FieldComponents::Spectral::SCALAR).fieldRange();
      for(fIt = fRange.first; fIt != fRange.second; ++fIt)
      {
      }

      // Make sure matrices are in compressed format
      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   void IBeta3DQGScalarEquation::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
   }
}
}
