/** \file BoussinesqShellVelocity.cpp
 *  \brief Source of the implementation of the Navier-Stokes equation in the Boussinesq spherical shell model
 */

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "Equations/Shell/Boussinesq/BoussinesqShellVelocity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralSelector.hpp"
#include "Equations/Tools/Vector2DEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqShellVelocity::BoussinesqShellVelocity(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
      // Vector equation has two components, ie Toroidal/Poloidal
      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqShellVelocity::~BoussinesqShellVelocity()
   {
   }

   void BoussinesqShellVelocity::initSpectralMatrices(const SharedSimulationBoundary spBcIds)
   {
      this->initSpectralMatrices2DEigen(spBcIds);
   }

   void BoussinesqShellVelocity::setCoupling()
   {
      // Get radial dimension
      int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get latitudinal dimension
      int nL = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get longitudinal dimension
      int nM = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

      //
      // Initialise  toroidal coupling information
      //
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::ONE,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::ONE);

      // General setup: prognostic equation, real solver, start from l = 1
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 1);

      // Set nonlinear flags: has NO nonlinear term, has NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::ONE);

      //
      // Initialise  poloidal coupling information
      //
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::TWO,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::TWO);

      // General setup: prognostic equation, real solver, start from l = 1
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 1);

      // Set nonlinear flags: has NO nonlinear term, has NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::TWO);

      //
      // Initialise common settings
      //

      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt = this->spectralRange();
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         // Set sizes of blocks and matrices
         ArrayI blockNs(nY);
         blockNs.setConstant(nX*nZ);
         ArrayI rhsCols(nY);
         rhsCols.setConstant(1);
         this->mCouplingInfos.find(*specIt)->second.setSizes(nY, blockNs, rhsCols);

         // Sort implicit fields
         this->mCouplingInfos.find(*specIt)->second.sortImplicitFields(eqId.first, *specIt);
      }
   }

   void BoussinesqShellVelocity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      throw Exception("Nonlinear term in vector equation not yet done");
   }

   void BoussinesqShellVelocity::setRequirements()
   {
      // Set velocity as equation unknown
      this->setName(PhysicalNames::VELOCITY);

      // Add velocity to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITY, FieldRequirement(false, true, false, false));
   }

   DecoupledZSparse BoussinesqShellVelocity::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return Vector2DEigenTools::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return Vector2DEigenTools::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return Vector2DEigenTools::boundaryRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void BoussinesqShellVelocity::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      throw Exception("setQuasiInverse for vector equation not yet done");
   }

   void BoussinesqShellVelocity::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k) const
   {
      throw Exception("setExplicitLinearBlock for vector equation not yet done");
   }

   void quasiInverseBlock(const BoussinesqShellVelocity& eq, SparseMatrix& mat)
   {
      throw Exception("Not yet implemented in vector equation");
   }

   void linearBlock(const BoussinesqShellVelocity& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
   }
      throw Exception("Not yet implemented in vector equation");

   void timeBlock(const BoussinesqShellVelocity& eq, DecoupledZSparse& mat, const MHDFloat k)
   {
      throw Exception("Not yet implemented in vector equation");
   }

   void boundaryBlock(const BoussinesqShellVelocity& eq, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat k)
   {
      throw Exception("Not yet implemented in vector equation");
   }

}
}
