/** 
 * @file BoussinesqBetaCylGVorticity.cpp
 * @brief Source of the implementation of the vorticity computation in the 3DQG beta model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/Beta3DQG/Boussinesq/BoussinesqBetaCylGVorticity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "SpectralOperators/Tools/SpectralBoxTools.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqBetaCylGVorticity::BoussinesqBetaCylGVorticity(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqBetaCylGVorticity::~BoussinesqBetaCylGVorticity()
   {
   }

   void BoussinesqBetaCylGVorticity::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: trivial equation, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to vorticity equation (self)
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Equation has explicit temperature
      infoIt.first->second.addExplicitField(PhysicalNames::STREAMFUNCTION,FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat = 0;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void BoussinesqBetaCylGVorticity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Doesn't require a nonlinear computation
      assert(false);
   }

   void BoussinesqBetaCylGVorticity::setRequirements()
   {
      // Set streamfunction as equation unknown
      this->setName(PhysicalNames::VORTICITY);

      // Set vorticity requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VORTICITY, FieldRequirement(true, true, true, true));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false));
   }

   void BoussinesqBetaCylGVorticity::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      // Quasi inverse is not required without nonlinear computations
      assert(false);
   }

   void BoussinesqBetaCylGVorticity::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      // Safety assert
      assert(compId == FieldComponents::Spectral::SCALAR);

      // PythonGenerator(mat, EXPLICIT, fieldId, compId, res, params, eigs, bcs)
   }

}
}
