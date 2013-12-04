/** 
 * @file BoussinesqShellVelocity.cpp
 * @brief Source of the implementation of the Navier-Stokes equation in the Boussinesq spherical shell model
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
#include "Equations/Shell/Boussinesq/BoussinesqShellVelocity.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"
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
         return EigenSelector::timeRow(*this, compId, matIdx);
      } else if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx);
      } else if(opId == IEquation::BOUNDARYROW)
      {
         return EigenSelector::boundaryRow(*this, compId, matIdx);
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

   void BoussinesqShellVelocity::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      //Safety assert
      assert(std::find(this->mSpectralIds.begin(), this->mSpectralIds.end(), compId) != this->mSpectralIds.end()); 

      linearBlock(*this, compId, mat, fieldId, eigs);
   }

   void quasiInverseBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, SparseMatrix& mat)
   {
      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nR);

      EigenSelector::KRSum blocks;

      // Set quasi-inverse for toroidal component
      if(compId == FieldComponents::Spectral::ONE)
      {
         blocks = spec1D.qDiff(2,0);

      // Set quasi-inverse for poloidal component
      } else if(compId == FieldComponents::Spectral::TWO)
      {
         blocks = spec1D.qDiff(2,0);

      } else
      {
         throw Exception("Unknown equation component ID for quasi-inverse operator!");
      }

      EigenSelector::computeKSum(mat, blocks);
   }

   void linearBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 2);
      MHDFloat l = eigs.at(0);
      MHDFloat m = eigs.at(1);

      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nR);

      EigenSelector::KZSum blocks;

      // Toroidal component : Velocity toroidal component
      if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         blocks.real() = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Toroidal component: Velocity poloidal component
      } else if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         blocks.real() = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Poloidal component : Velocity toroidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         blocks.real() = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Poloidal component: Velocity poloidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         blocks.real() = Spectral::SphericalHarmonicOperator::qLaplacian(spec1D, l, m, 2);

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID and component ID combination for linear operator!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs);
   }

   void timeBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs)
   {
      assert(eigs.size() == 2);
      MHDFloat l = eigs.at(0);
      MHDFloat m = eigs.at(1);

      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Get R
      int nR = eq.unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);

      // Create spectral operator
      Spectral::OperatorSelector<Dimensions::Simulation::SIM1D>::Type spec1D(nR);

      EigenSelector::KZSum blocks;

      // Time operator for toroidal component
      if(compId == FieldComponents::Spectral::ONE)
      {
         blocks.real() = spec1D.qDiff(2,0);

      // Time operator for poloidal component
      } else if(compId == FieldComponents::Spectral::TWO)
      {
         blocks.real() = spec1D.qDiff(2,0);

      } else
      {
         throw Exception("Unknown field ID and component ID combination for time operator!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs);
   }

   void boundaryBlock(const BoussinesqShellVelocity& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 2);
      MHDFloat l = eigs.at(0);
      MHDFloat m = eigs.at(1);

      //Safety assert
      assert(std::find(eq.spectralRange().first, eq.spectralRange().second, compId) != eq.spectralRange().second); 

      // Toroidal component: Boundary condition for the velocity toroidal component
      if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Toroidal component: Boundary condition for the velocity poloidal component
      } else if(compId == FieldComponents::Spectral::ONE && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Poloidal component: Boundary condition for the velocity toroidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::ONE)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Poloidal component: Boundary condition for the velocity poloidal component
      } else if(compId == FieldComponents::Spectral::TWO && fieldId.first == PhysicalNames::VELOCITY && fieldId.second == FieldComponents::Spectral::TWO)
      {
         coeffs.push_back(1.0);
         bcIdx.push_back(Boundary::BCIndex(Boundary::INDEPENDENT));

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID and component ID combination for boundary condition operator!");
      }
   }

}
}
