/** 
 * @file TestTFTExactForwardScalar.cpp
 * @brief Source of the implementation of a test equation for the TFT scheme with exact known scalar physical space solution
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
#include "Equations/Tests/TestTFTExactForwardScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFTExactForwardScalar::TestTFTExactForwardScalar()
      : IScalarEquation(SharedEquationParameter(new EquationParameters()))
   {
   }

   TestTFTExactForwardScalar::~TestTFTExactForwardScalar()
   {
   }

   void TestTFTExactForwardScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestTFTExactForwardScalar::setCoupling()
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

      // General setup: trivial equation, real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, true, 0);

      // Set nonlinear flags: HAS nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(true, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Standalone equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Set sizes of blocks and matrices
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void TestTFTExactForwardScalar::setRequirements()
   {
      // Set simple requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   void TestTFTExactForwardScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Doesn't require a nonlinear computation
      assert(false);
   }

}
}
