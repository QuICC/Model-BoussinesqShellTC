/** \file BoussinesqBetaCylGVertical.cpp
 *  \brief Source of the implementation of the vertical velocity equation in the 3DQG beta model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVertical.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGVertical::BoussinesqBetaCylGVertical(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGVertical::~BoussinesqBetaCylGVertical()
   {
   }

   void BoussinesqBetaCylGVertical::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices1DPeriodic(spBcIds);
   }

   void BoussinesqBetaCylGVertical::setCoupling()
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

      // General setup: prognostic equation, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 0);

      // set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(true, true);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // 
      //  WARNING: the order is important as it determines the field index!
      //

      // Equation is coupled to streamfunction equation
      infoIt.first->second.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR, false);
      // Equation is coupled to vertical velocity equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR, true);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 
   }

   void BoussinesqBetaCylGVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 0.0);
   }

   void BoussinesqBetaCylGVertical::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   DecoupledZSparse BoussinesqBetaCylGVertical::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return timeRow1DPeriodic(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
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

   void BoussinesqBetaCylGVertical::setQuasiInverse(SparseMatrix& mat) const
   {
      quasiInverseBlock(*this, mat);
   }

   void BoussinesqBetaCylGVertical::setExplicitLinearBlock(DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      linearBlock(*this, mat, fieldId, k);
   }

   void quasiInverseBlock(const BoussinesqBetaCylGVertical& eq, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nX = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nX);
      Spectral::SpectralSelector<Dimensions::Simulation::SIM3D>::OpType spec3D(nZ);

      /// - Vertical velocity equation: \f$ \left(D_x^{-2} \otimes D_Z^{-1}\right) \f$
      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), mat);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqBetaCylGVertical& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
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

      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);

      /// - Streamfunction : \f$ -\frac{1}{\Gamma^2}\left(D_x^{-2} \otimes I_Z^{-1}\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         SparseMatrix tmp = -(1./(Gamma*Gamma))*spec3D.id(1);
         Eigen::kroneckerProduct(tmp, spec1D.qDiff(2,0), mat.first);

      /// - Vertical velocity : \f$ \left(D_x^{-2}\nabla_\perp^{2}\otimes D_Z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         Eigen::kroneckerProduct(spec3D.qDiff(1,0), Spectral::PeriodicOperator::qLaplacian2D(spec1D, k_, 2), mat.first);

      /// - Temperature : \f$ \left(0_x \otimes 0_Z\right) \f$
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         //
         // Nothing to do for an empty sparse matrix
         //

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }
   }

   void timeBlock(const BoussinesqBetaCylGVertical& eq, DecoupledZSparse& mat, const MHDFloat k)
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

      // Set time matrices (kronecker(A,B,out) => out = A(i,j)*B)
      Eigen::kroneckerProduct(spec3D.qDiff(1,0), spec1D.qDiff(2,0), mat.first);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const BoussinesqBetaCylGVertical& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Rescale wave number to [-1, 1]
      MHDFloat k_ = k/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = eq.eqParams().nd(NonDimensional::CHI);

      /// <b>Boundary operators for the vertical velocity equation</b>: \f$\left(BC_x \otimes I_Z\right)\f$ and \f$\left(D_x^{-2} \otimes BC_Z\right)\f$
      int pX = 2;
      int pZ = 0;

      // Set boundary condition prefactors
      MHDFloat cX = 1.0;
      MHDFloat cZ = k_*std::tan((MathConstants::PI/180.)*chi)/Gamma;

      // Compute boundary block operator
      boundaryBlock1DPeriodic(eq, mat, fieldId, pX, pZ, cX, cZ);
   }

}
}
