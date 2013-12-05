/** 
 * @file TestFFFDiffusion3D.cpp
 * @brief Source of the implementation of the FFF test equation for 3D diffusion
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
#include "Equations/Tests/TestFFFDiffusion3D.hpp"

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

   TestFFFDiffusion3D::TestFFFDiffusion3D(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
   }

   TestFFFDiffusion3D::~TestFFFDiffusion3D()
   {
   }

   void TestFFFDiffusion3D::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestFFFDiffusion3D::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: prognostic equation real solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::PROGNOSTIC, false, 0);

      // Set nonlinear flags: has nonlinear term, has quasi-inverse
      infoIt.first->second.setNonlinear(false, false);

      // Set source flags: NO source term
      infoIt.first->second.setSource(false);

      // Equation is coupled to itself
      infoIt.first->second.addImplicitField(eqId.first, FieldComponents::Spectral::SCALAR);

      // Set mininal matrix coupling
      int nMat = 0;
      ArrayI blockNs;
      ArrayI rhsCols;
      EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
      infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void TestFFFDiffusion3D::setRequirements()
   {
      // Add requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   DecoupledZSparse TestFFFDiffusion3D::operatorRow(const IEquation::OperatorRowId opId, FieldComponents::Spectral::Id compId, const int matIdx, const bool hasBoundary) const
   {
      if(opId == IEquation::TIMEROW)
      { 
         return EigenSelector::timeRow(*this, compId, matIdx, hasBoundary);
      } else if(opId == IEquation::LINEARROW)
      {
         return EigenSelector::linearRow(*this, compId, matIdx, hasBoundary);
      } else
      {
         throw Exception("Unknown operator row ID");
      }
   }

   void linearBlock(const TestFFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 3);
      // Rescale wave number to [-1, 1]
      MHDFloat k_ = eigs.at(0)/2.0;
      MHDFloat l_ = eigs.at(1)/2.0;
      MHDFloat m_ = eigs.at(2)/2.0;

      EigenSelector::KZSum blocks;

      // Setup 3D diffusion
      if(fieldId.first == eq.name())
      {
         blocks.real() = Spectral::BoxTools::qLaplacian3D(k_, l_, m_);

      // Unknown field
      } else
      {
         throw Exception("Unknown field ID for linear operator!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void timeBlock(const TestFFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, const bool hasBoundary)
   {
      assert(eigs.size() == 3);

      EigenSelector::KZSum blocks;

      if(fieldId.first == eq.name())
      {
         blocks.real() = 1.0;
      } else
      {
         throw Exception("Multiple field in time integration not implemented yet!");
      }

      EigenSelector::constrainBlock(eq, compId, mat, fieldId, blocks, eigs, hasBoundary);
   }

   void boundaryBlock(const TestFFFDiffusion3D& eq, FieldComponents::Spectral::Id compId, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs, std::vector<MHDFloat>& coeffs, std::vector<Boundary::BCIndex>& bcIdx)
   {
      assert(eigs.size() == 3);

      throw Exception("There are not boundary conditions to impose!");
   }

}
}
