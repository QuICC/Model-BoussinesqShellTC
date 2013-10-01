/** 
 * @file BoussinesqShellTransport.cpp
 * @brief Source of the implementation of the transport equation in the spherical shell model
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
#include "Equations/Shell/Boussinesq/BoussinesqShellTransport.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "TypeSelectors/EquationToolsSelector.hpp"
#include "SpectralOperators/SphericalHarmonicOperator.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqShellTransport::BoussinesqShellTransport(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqShellTransport::~BoussinesqShellTransport()
   {
   }

   void BoussinesqShellTransport::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: prognostic equation, real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 0);

      // Set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Set index type: use MODE
      infoIt.first->second.setIndexType(CouplingInformation::MODE);

      // Equation is coupled to temperature equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);

      //
      // Initialise common dimensional settings
      //

      // Set mininal matrix coupling
      int nMat;
      ArrayI blockNs;
      ArrayI rhsCols;
      EquationToolsType::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 
   }

   void BoussinesqShellTransport::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      throw Exception("Nonlinear term in spherical shell transport equation not done yet");
   }

   void BoussinesqShellTransport::setRequirements()
   {
      // Set temperatur as equation unknown
      this->setName(PhysicalNames::TEMPERATURE);

      // Add temperature to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse BoussinesqShellTransport::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return EquationToolsType::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return EquationToolsType::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return EquationToolsType::boundaryRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void BoussinesqShellTransport::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqShellTransport::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      linearBlock(*this, compId, mat, fieldId, k);
   }

   void quasiInverseBlock(const BoussinesqShellTransport& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Set quasi-inverse
      mat = spec1D.qDiff(2,0);

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqShellTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Initialise output matrices
      mat.first.resize(nR,nR);
      mat.second.resize(nR,nR);

      // Temperature
      if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Velocity toroidal component
      } else if(fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Velocity poloidal component
      } else if(fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID and component ID combination for linear operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void timeBlock(const BoussinesqShellTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Initialise output matrices
      mat.first.resize(nR,nR);
      mat.second.resize(nR,nR);

      // Time operator
      mat.first = spec1D.qDiff(2,0);

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);

   }

   void boundaryBlock(const BoussinesqShellTransport& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR); 

      // Storage for the boundary condition constant factor
      MHDFloat cR;

      // Boundary condition for the temperature
      if(fieldId.first == PhysicalNames::TEMPERATURE)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Boundary condition for the velocity toroidal component
      } else if(fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Boundary condition for the velocity poloidal component
      } else if(fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID and component ID combination for boundary condition operator!");
      }

      // Compute boundary block operator
      EquationToolsType::boundaryBlock(eq, FieldComponents::Spectral::SCALAR, mat, fieldId, cR);
   }

}
}
