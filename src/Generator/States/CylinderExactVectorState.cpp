/** 
 * @file CylinderExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a cylinder
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
#include "Generator/States/CylinderExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   CylinderExactVectorState::CylinderExactVectorState(const std::string& pyName, SharedEquationParameters spEqParams)
      : IVectorEquation(pyName,spEqParams), mTypeId(CONSTANT)
   {
   }

   CylinderExactVectorState::~CylinderExactVectorState()
   {
   }

   void CylinderExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      this->mSpectralIds.push_back(FieldComponents::Spectral::ONE);
      this->mSpectralIds.push_back(FieldComponents::Spectral::TWO);

      // Set the variable requirements
      this->setRequirements();
   }

   void CylinderExactVectorState::setStateType(const CylinderExactVectorState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void CylinderExactVectorState::setCoupling()
   {
      SpectralComponent_range specRange = this->spectralRange();
      for(SpectralComponent_iterator it = specRange.first; it != specRange.second; ++it)
      {
         // Initialise coupling information
         std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
         infoIt = this->mCouplingInfos.insert(std::make_pair(*it,CouplingInformation()));
         SpectralFieldId eqId = std::make_pair(this->name(), *it);

         // General setup: first complex solver, complex solver, start from m = 0
         infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

         // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
         infoIt.first->second.setNonlinear(true, false);

         // Set source flags: has source term
         infoIt.first->second.setSource(false);

         // Equation is coupled to itself
         infoIt.first->second.addImplicitField(eqId.first, *it);

         // Set mininal matrix coupling
         int nMat = 0;
         ArrayI blockNs;
         ArrayI rhsCols;
         EigenSelector::makeMinimalCoupling(this->unknown().dom(0).spRes(), nMat, blockNs, rhsCols);
         infoIt.first->second.setSizes(nMat, blockNs, rhsCols); 

         // Sort implicit fields
         infoIt.first->second.sortImplicitFields(eqId.first, *it);
      }
   }

   void CylinderExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == CONSTANT)
      {
         int nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array zGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               for(int iZ = 0; iZ < nZ; ++iZ)
               {
                  rNLComp.setPoint(1.0, iZ,iTh,iR);
               }
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex CylinderExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return MHDComplex(0,0);
    }

   void CylinderExactVectorState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

   void CylinderExactVectorState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // PythonGenerator(mat, QI, fieldId, compId, res, params, eigs, bcs)
   }

   void CylinderExactVectorState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      // PythonGenerator(mat, EXPLICIT, fieldId, compId, res, params, eigs, bcs)
   }

}
}
