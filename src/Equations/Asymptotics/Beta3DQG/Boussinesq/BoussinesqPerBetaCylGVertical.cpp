/** 
 * @file BoussinesqPerBetaCylGVertical.cpp
 * @brief Source of the implementation of the vertical velocity equation in the 3DQG beta model with periodic radius
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqPerBetaCylGVertical.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "PhysicalOperators/StreamAdvection.hpp"
#include "SpectralOperators/PeriodicOperator.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqPerBetaCylGVertical::BoussinesqPerBetaCylGVertical(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqPerBetaCylGVertical::~BoussinesqPerBetaCylGVertical()
   {
   }

   void BoussinesqPerBetaCylGVertical::setCoupling()
   {
      // Set coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: prognostic equation, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 1);

      // set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to streamfunction equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Equation is coupled to vertical velocity equation
      infoIt.first->second.addImplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EquationToolsType::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqPerBetaCylGVertical::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      /// 
      /// Computation of the jacobian:
      ///   \f$ \left(\nabla^{\perp}\psi\cdot\nabla_{\perp}\right)w\f$
      ///
      Physical::StreamAdvection::set(rNLComp, this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad(), this->unknown().dom(0).grad(), 0.0);
   }

   void BoussinesqPerBetaCylGVertical::setRequirements()
   {
      // Set vertical velocity as equation unknown
      this->setName(PhysicalNames::VELOCITYZ);

      // Add vertical velocity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYZ, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, false, false, true));
   }

   DecoupledZSparse BoussinesqPerBetaCylGVertical::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return timeRow2DPeriodic(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return linearRow2DPeriodic(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return boundaryRow2DPeriodic(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void BoussinesqPerBetaCylGVertical::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqPerBetaCylGVertical::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      linearBlock(*this, compId, mat, fieldId, k);
   }

   void quasiInverseBlock(const BoussinesqPerBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      // Get X and Z dimensions
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Set quasi-inverse operator of streamfunction equation multiplication matrix (kronecker(A,B,out) => out = A(i,j)*B)
      mat = spec1D.qDiff(1,0);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqPerBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Get X and Z dimensions
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Initialise output matrices
      mat.real().resize(nZ,nZ);
      mat.imag().resize(nZ,nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat kX_ = kX/2.;
      MHDFloat kY_ = kY/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Ra = eq.eqParams().nd(NonDimensional::RAYLEIGH);
      MHDFloat Pr = eq.eqParams().nd(NonDimensional::PRANDTL);

      /// - Streamfunction : \f$ \left(D_x^{-4}\nabla_\perp^{4} \otimes D_Z^{-1}\right) \f$
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.real() = -(1./(Gamma*Gamma))*spec1D.id(1);

      /// - Vertical velocity : \f$ \left(D_x^{-4} \otimes I_z^{-1}\right) \f$
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Build linear operator (kronecker(A,B,out) => out = A(i,j)*B)
         mat.real() = Spectral::PeriodicOperator::laplacian2D(kX_, kY_)*spec1D.qDiff(1,0);

      /// - Temperature : \f$ i \frac{k}{2}\frac{1}{16}\frac{Ra}{Pr}\left( D_x^{-4} \otimes D_Z^{-1}\right) \f$
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

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void timeBlock(const BoussinesqPerBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat k)
   {
      // Get X and Z dimensions
      int nZ = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operators
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nZ);

      // Initialise output matrices
      mat.real().resize(nZ,nZ);
      mat.imag().resize(nZ,nZ);

      // Rescale wave number to [-1, 1]
      MHDFloat kX_ = kX/2.;
      MHDFloat kY_ = kY/2.;

      // Set time matrix (kronecker(A,B,out) => out = A(i,j)*B)
      mat.real() = spec1D.qDiff(1,0);

      // Prune matrices for safety
      mat.real().prune(1e-32);
      mat.imag().prune(1e-32);
   }

   void boundaryBlock(const BoussinesqPerBetaCylGVertical& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      // Rescale wave number to [-1, 1]
      MHDFloat kY_ = kY/2.;

      // Get Physical parameters Ra, Pr, Gamma, chi
      MHDFloat Gamma = eq.eqParams().nd(NonDimensional::GAMMA);
      MHDFloat chi = eq.eqParams().nd(NonDimensional::CHI);

      int pZ = 0;

      MHDFloat cX;
      MHDFloat cZ;

      // Boundary condition for the streamfunction
      if(fieldId.first == PhysicalNames::STREAMFUNCTION)
      {
         // Set boundary condition prefactors
         cZ = k_*std::tan((MathConstants::PI/180.)*chi)/Gamma;

      // Boundary condition for the vertical velocity
      } else if(fieldId.first == PhysicalNames::VELOCITYZ)
      {
         // Set boundary condition prefactors
         cZ = 1.0;

      // Boundary condition for the temperature
      } else if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Set boundary condition prefactors
         cZ = 0.0;

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for boundary operator!");
      }

      // Compute boundary block operator
      boundaryBlock2DPeriodic(eq, mat, fieldId, pZ, cZ);
   }

}
}
