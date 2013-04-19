/** \file Beta3DQGStreamfunction.cpp
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
#include "Equations/Asymptotics/Beta3DQG/Beta3DQGStreamfunction.hpp"

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

   Beta3DQGStreamfunction::Beta3DQGStreamfunction(SharedIEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Beta3DQGStreamfunction::~Beta3DQGStreamfunction()
   {
   }

   void Beta3DQGStreamfunction::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp) const
   {
      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)\nabla^2_{\perp}\psi\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->unknown().dom(0).grad(), this->scalar(PhysicalNames::VORTICITYZ).dom(0).grad(), 1.0);
   }

   void Beta3DQGStreamfunction::computeLinear(Datatypes::SpectralScalarType& rRHS) const
   {  
      ///
      /// Compute \f$-\frac{1}{16}\frac{Ra}{Pr}\partial_y\overline{T} = -\frac{1}{16}\frac{Ra}{Pr} i m/2 \overline{T}\f$
      ///

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Compute Ra/(2 16 Pr)
      MHDComplex c = -this->eqParams().nd(NonDimensional::RAYLEIGH)/(32.*this->eqParams().nd(NonDimensional::PRANDTL))*boxScale*MathConstants::cI;

      // Get size of dealiased output (at this stage data has still dealiasing rows)
      int dealiasedRows = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Loop over m
      MHDFloat m_;
      int nSlice = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      for(int m = 0; m < nSlice; m++)
      {
         m_ = static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(m));

         rRHS.addSlice((m_*c)*this->scalar(PhysicalNames::TEMPERATURE).dom(0).perturbation().slice(m), m, dealiasedRows);
      }
   }

   void Beta3DQGStreamfunction::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::STREAMFUNCTION);

      // Set streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, true, true));

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITYZ, FieldRequirement(true, true, false, true));

      // Add temperature requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, false, false, true));
   }

   void Beta3DQGStreamfunction::setCoupling()
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

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const DecoupledZMatrix& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Get right wave number
      MHDFloat m_ = boxScale*0.5*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }

   void Beta3DQGStreamfunction::timestepOutput(FieldComponents::Spectral::Id id, const MatrixZ& storage, const int matIdx, const int start)
   {
      // Call basic implementation
      IScalarEquation::timestepOutput(id, storage, matIdx, start);

      // Get the box scale
      MHDFloat boxScale = this->unknown().dom(0).spRes()->sim()->boxScale(Dimensions::Simulation::SIM2D);

      // Get right wave number
      MHDFloat m_ = boxScale*0.5*static_cast<MHDFloat>(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(matIdx));

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL));

      ///
      /// Compute the vertical vorticity: \f$\zeta = \nabla^2\psi\f$
      ///
      this->rScalar(PhysicalNames::VORTICITYZ).rDom(0).rPerturbation().setSlice(Spectral::PeriodicOperator::laplacian2D(spec1D, m_, 0)*this->unknown().dom(0).perturbation().slice(matIdx), matIdx);
   }

   DecoupledZSparse Beta3DQGStreamfunction::linearRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

   DecoupledZSparse Beta3DQGStreamfunction::timeRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

   DecoupledZSparse Beta3DQGStreamfunction::boundaryRow(FieldComponents::Spectral::Id comp, const int matIdx) const
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

   void Beta3DQGStreamfunction::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
   }
}
}
