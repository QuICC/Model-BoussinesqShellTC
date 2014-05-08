/** 
 * @file SphereExactScalarState.cpp
 * @brief Source of the implementation of the equation to generate an exact scalar solution in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//
#include <tr1/cmath>

// External includes
//

// Class include
//
#include "Generator/States/SphereExactScalarState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/SpectralOperatorSelector.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SphereExactScalarState::SphereExactScalarState(const std::string& pyName, SharedEquationParameters spEqParams)
      : IScalarEquation(pyName,spEqParams), mTypeId(CONSTANT)
   {
   }

   SphereExactScalarState::~SphereExactScalarState()
   {
   }

   void SphereExactScalarState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void SphereExactScalarState::setStateType(const SphereExactScalarState::StateTypeId id)
   {
      this->mTypeId = id;
   }

   void SphereExactScalarState::setHarmonicOptions(const std::vector<std::tr1::tuple<int,int,MHDComplex> >& modes)
   {
      this->mSHModes = modes;
   }

   void SphereExactScalarState::setCoupling()
   {
      // Initialise coupling information
      std::pair<std::map<FieldComponents::Spectral::Id, CouplingInformation>::iterator,bool> infoIt;
      infoIt = this->mCouplingInfos.insert(std::make_pair(FieldComponents::Spectral::SCALAR,CouplingInformation()));
      SpectralFieldId eqId = std::make_pair(this->name(), FieldComponents::Spectral::SCALAR);

      // General setup: first complex solver, complex solver, start from m = 0
      infoIt.first->second.setGeneral(CouplingInformation::TRIVIAL, false, 0);

      // Set nonlinear flags: NO nonlinear term, NO quasi-inverse
      infoIt.first->second.setNonlinear(true, false);

      // Set source flags: has source term
      infoIt.first->second.setSource(false);

      // Set index type: use MODE
      infoIt.first->second.setIndexType(CouplingInformation::MODE);

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

   void SphereExactScalarState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Physical::SCALAR);

      if(this->mTypeId == CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(this->mTypeId == HARMONIC)
      {
         int nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         Array funcPh(nPh);
         MHDFloat funcR(nR);
         MHDFloat funcTh;
         std::vector<std::tr1::tuple<int,int,MHDComplex> >::const_iterator it;
         rNLComp.rData().setConstant(0);
         for(int iR = 0; iR < nR; ++iR)
         {
            for(int iTh = 0; iTh < nTh; ++iTh)
            {
               for(it = this->mSHModes.begin(); it != this->mSHModes.end(); ++it)
               {
                  int l = std::tr1::get<0>(*it);
                  int m = std::tr1::get<1>(*it);
                  MHDFloat re = std::tr1::get<2>(*it).real();
                  MHDFloat im = std::tr1::get<2>(*it).imag();

                  funcR = std::pow(rGrid(iR),l);

                  // Spherical harmonic Y_l^m
                  funcPh = re*(static_cast<MHDFloat>(m)*phGrid).array().cos() + im*(static_cast<MHDFloat>(m)*phGrid).array().sin();
                  funcTh = std::tr1::sph_legendre(l,m, thGrid(iTh));
                  rNLComp.addProfile(funcPh*funcR*funcTh,iTh,iR);
               }
            }
         }

      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    MHDComplex SphereExactScalarState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      return MHDComplex(0,0);
    }

   void SphereExactScalarState::setRequirements()
   {
      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, false, false));
   }

   void SphereExactScalarState::setQuasiInverse(FieldComponents::Spectral::Id compId, SparseMatrix& mat) const
   {
      // PythonGenerator(mat, QI, fieldId, compId, res, params, eigs, bcs)
   }

   void SphereExactScalarState::setExplicitLinearBlock(FieldComponents::Spectral::Id compId, DecoupledZSparse& mat, const SpectralFieldId fieldId, const std::vector<MHDFloat>& eigs) const
   {
      // PythonGenerator(mat, EXPLICIT, fieldId, compId, res, params, eigs, bcs)
   }

}
}
