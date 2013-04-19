/** \file Beta3DQGVertical.cpp
 *  \brief Source of the implementation of the vertical velocity equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//
#include <Eigen/Sparse>
#include <Eigen/KroneckerProduct>

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGVertical.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGSystem.hpp"

#include <iostream>

namespace GeoMHDiSCC {

namespace Equations {

   Beta3DQGVertical::Beta3DQGVertical(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGVertical::~Beta3DQGVertical()
   {
   }

   void Beta3DQGVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 1.0);
   }

   void Beta3DQGVertical::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGVertical::setCoupling()
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

   DecoupledZSparse Beta3DQGVertical::linearRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   DecoupledZSparse Beta3DQGVertical::timeRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

      std::cerr << matrixRow.first << std::endl;
      std::cerr << matrixRow.second << std::endl;

      return matrixRow;
   }

   DecoupledZSparse Beta3DQGVertical::boundaryRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

      matrixRow.first.makeCompressed();
      matrixRow.second.makeCompressed();

      return matrixRow;
   }

   void Beta3DQGVertical::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
   }
}
}
