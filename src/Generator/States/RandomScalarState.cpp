/** \file RandomScalarState.cpp
 *  \brief Source of the implementation of the general random scalar state equation
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Generator/States/RandomScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   RandomScalarState::RandomScalarState(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   RandomScalarState::~RandomScalarState()
   {
   }

   void RandomScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void RandomScalarState::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DPeriodic(spBcIds);
   }

   void RandomScalarState::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first complex solver, complex solver, start from m = 0
      infoIt.first->second.setGeneral(0, false, 0);

      // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: has source term
      infoIt.first->second.setSource(false);

      // 
      //  WARNING: the order is important as it determines the field index!
      //

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 
   }

   void RandomScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);
   }

   void RandomScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse RandomScalarState::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::LINEARROW)
      {
         return linearRow1DPeriodic(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return boundaryRow1DPeriodic(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void RandomScalarState::setQuasiInverse(SparseMatrix& mat) const
   {
      Equations::quasiInverseBlock(*this, mat);
   }

   void RandomScalarState::setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      Equations::linearBlock(*this, mat, fieldId, k);
   }

   void quasiInverseBlock(const RandomScalarState& eq, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0), mat);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const RandomScalarState& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      // Initialise output matrices
      mat.first.resize(nX*nZ,nX*nZ);
      mat.second.resize(nX*nZ,nX*nZ);

      // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.id(0), spec1D.id(0), mat.first);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const RandomScalarState& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      int pX = 0;
      int pZ = 0;

      // Set boundary condition prefactors
      MHDFloat cX = 1.0;
      MHDFloat cZ = 1.0;

      // Compute boundary block operator
      boundaryBlock1DPeriodic(eq, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
