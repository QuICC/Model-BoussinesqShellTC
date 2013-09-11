/** 
 * @file BoussinesqShellVelocity.cpp
 * @brief Source of the implementation of the Navier-Stokes equation in the Boussinesq spherical shell model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
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
#include "Equations/Tools/Equation2DEigenTools.hpp"
#include "SpectralOperators/SphericalHarmonicOperator.hpp"

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
      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator compIt;
      for(compIt = specRange.first; compIt != specRange.second; ++compIt)
      {
         this->initSpectralMatrices2DEigen(spBcIds, *compIt);
      }
   }

   void BoussinesqShellVelocity::setCoupling()
   {
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

      // Set index type: use MODE
      infoIt.first->second.setIndexType(CouplingInformation::MODE);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::ONE);

      //
      // Initialise  poloidal coupling information
      //
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::TWO,CouplingInformation()));
      eqId = std::make_pair(this->name(), FieldComponents::Spectral::TWO);

      // General setup: prognostic equation, real solver, start from l = 1
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, true, 1);

      // Set nonlinear flags: has NO nonlinear term, has NO quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Set index type: use MODE
      infoIt.first->second.setIndexType(CouplingInformation::MODE);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::TWO);

      //
      // Initialise common dimensional settings
      //
      
      // Get radial dimension
      int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get latitudinal dimension
      int nL = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

      SpectralComponent_range specRange = this->spectralRange();
      SpectralComponent_iterator specIt;
      for(specIt = specRange.first; specIt != specRange.second; ++specIt)
      {
         // Set sizes of blocks and matrices
         ArrayI blockNs(nL);
         blockNs.setConstant(nR);
         ArrayI rhsCols(nL);
         for(int i = 0; i < rhsCols.size(); i++)
         {
            rhsCols(i) = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(i);
         }
         this->mCouplingInfos.find(*specIt)->second.setSizes(nL, blockNs, rhsCols);

         // Sort implicit fields
         this->mCouplingInfos.find(*specIt)->second.sortImplicitFields(eqId.first, *specIt);
      }
   }

   void BoussinesqShellVelocity::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      throw Exception("Nonlinear term in spherical shell toroidal/poloidal Navier-Stokes equation not done yet");
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
         return Equation2DEigenTools::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return Equation2DEigenTools::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return Equation2DEigenTools::boundaryRow(*this, compId, matIdx);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void BoussinesqShellVelocity::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      //Safety assert
      assert(std::find(this->mSpectralIds.begin(), this->mSpectralIds.end(), compId) != this->mSpectralIds.end()); 

      quasiInverseBlock(*this, compId, mat);
   }

   void BoussinesqShellVelocity::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat l, const MHDFloat m) const
   {
      //Safety assert
      assert(std::find(this->mSpectralIds.begin(), this->mSpectralIds.end(), compId) != this->mSpectralIds.end()); 

      linearBlock(*this, compId, mat, fieldId, l, m);
   }

   void quasiInverseBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Set quasi-inverse for toroidal component
      if(compId == FieldComponents::Spectral::ONE)
      {
         mat = spec1D.qDiff(2,0);

      // Set quasi-inverse for poloidal component
      } else if(compId == FieldComponents::Spectral::TWO)
      {
         mat = spec1D.qDiff(2,0);

      } else
      {
         throw Exception("Unknown equation component ID for quasi-inverse operator!");
      }

      // Prune matrices for safety
      mat.prune(1e-32);
   }

   void linearBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Initialise output matrices
      mat.first.resize(nR,nR);
      mat.second.resize(nR,nR);

      // Toroidal component : Velocity toroidal component
      if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Toroidal component: Velocity poloidal component
      } else if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Poloidal component : Velocity toroidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         mat.first = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Poloidal component: Velocity poloidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
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

   void timeBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::SpectralSelector<Dimensions::Simulation::SIM1D>::OpType spec1D(nR);

      // Initialise output matrices
      mat.first.resize(nR,nR);
      mat.second.resize(nR,nR);

      // Time operator for toroidal component
      if(compId == FieldComponents::Spectral::ONE)
      {
         mat.first = spec1D.qDiff(2,0);

      // Time operator for poloidal component
      } else if(compId == FieldComponents::Spectral::TWO)
      {
         mat.first = spec1D.qDiff(2,0);

      } else
      {
         throw Exception("Unknown field ID and component ID combination for time operator!");
      }

      // Prune matrices for safety
      mat.first.prune(1e-32);
      mat.second.prune(1e-32);
   }

   void boundaryBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const MHDFloat l, const MHDFloat m)
   {
      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Storage for the boundary condition constant factor
      MHDFloat cR;

      // Toroidal component: Boundary condition for the velocity toroidal component
      if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Toroidal component: Boundary condition for the velocity poloidal component
      } else if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Poloidal component: Boundary condition for the velocity toroidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Poloidal component: Boundary condition for the velocity poloidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         // Set boundary condition prefactors
         cR = 1.0;

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID and component ID combination for boundary condition operator!");
      }

      // Compute boundary block operator
      Equation2DEigenTools::boundaryBlock2DEigen(eq, compId, mat, fieldId, cR);
   }

}
}
