/** 
 * @file TiltedScalarFieldVisualizer.cpp
 * @brief Source of the implementation of the tilted scalar field visualizer
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
#include "Generator/Visualizers/TiltedScalarFieldVisualizer.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TiltedScalarFieldVisualizer::TiltedScalarFieldVisualizer(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams), mViewField(true), mViewGradient(false)
   {
   }

   TiltedScalarFieldVisualizer::~TiltedScalarFieldVisualizer()
   {
   }

   void TiltedScalarFieldVisualizer::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TiltedScalarFieldVisualizer::setDataField(const PhysicalNames::Id name)
   {
      this->mDataField = name;
   }

   void TiltedScalarFieldVisualizer::setFields(const bool viewField, const bool viewGradient)
   {
      this->mViewField = viewField;

      this->mViewGradient = viewGradient;
   }

   void TiltedScalarFieldVisualizer::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::WRAPPER, 0, true, false);
   }

   void TiltedScalarFieldVisualizer::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {  
      // Get paramters
      MHDFloat eta3 = std::cos((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));
      MHDFloat eta2 = std::sin((Math::PI/180.)*this->eqParams().nd(NonDimensional::THETA));

      // Initialize FFT
      Transform::Fft::FftSelector   transform;
      Transform::Fft::FftSelector::SharedSetupType spSetup = std::tr1::static_pointer_cast<Transform::Fft::FftSelector::SetupType>(this->unknown().dom(0).spRes()->spTransformSetup(Dimensions::Transform::TRA3D));
      transform.init(spSetup);

      // Compute forward transform
      MatrixZ tmp(spSetup->bwdSize(), spSetup->howmany());
      transform.integrate(tmp, this->scalar(this->mDataField).dom(0).phys().data(), Transform::Fft::FftSelector::IntegratorType::INTG);

      // Get Z grid
      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ;
      int nI;
      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);

      MHDFloat k_;
      int m = 0;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = (1.0 - gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK)))/2.0;
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            nI = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DATB1D>(iK);
            for(int iI = 0; iI < nI; ++iI)
            {
               tmp(iI, m) = std::exp(MHDComplex(0.0, -k_*iI*(eta2/eta3)))*tmp(iI,m);
            }
            m++;
         }
      }

      // Compute backward transform
      transform.project(rNLComp.rData(), tmp, Transform::FftIds::Projectors::PROJ);
   }

   void TiltedScalarFieldVisualizer::useNonlinear(const Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId)
   {  
      this->rUnknown().rDom(0).rPhys().rData() = rNLComp.data();
   }

   void TiltedScalarFieldVisualizer::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, this->mViewField, this->mViewGradient));

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->mDataField, FieldRequirement(true, true, this->mViewField, this->mViewGradient));
   }

}
}
