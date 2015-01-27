/** 
 * @file ShellExactVectorState.cpp
 * @brief Source of the implementation of the equation to generate an exact vector solution in a shell
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
#include "Generator/States/ShellExactVectorState.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   ShellExactVectorState::ShellExactVectorState(SharedEquationParameters spEqParams)
      : IVectorEquation(spEqParams)
   {
   }

   ShellExactVectorState::~ShellExactVectorState()
   {
   }

   void ShellExactVectorState::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void ShellExactVectorState::setStateType(const ShellExactStateIds::Id id)
   {
      this->mTypeId = id;
   }

   void ShellExactVectorState::setHarmonicOptions(const FieldComponents::Spectral::Id compId, const std::vector<ShellExactVectorState::HarmonicModeType>& modes)
   {
      this->mSHModes.insert(std::make_pair(compId, modes));
   }

   void ShellExactVectorState::setCoupling()
   {
      if(FieldComponents::Spectral::ONE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::ONE, CouplingInformation::TRIVIAL, 0, true, false, false, false);
         this->setExplicitTiming(FieldComponents::Spectral::ONE, ExplicitTiming::LINEAR);
      }

      if(FieldComponents::Spectral::TWO != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::TWO, CouplingInformation::TRIVIAL, 0, true, false, false, false);
         this->setExplicitTiming(FieldComponents::Spectral::TWO, ExplicitTiming::LINEAR);
      }

      if(FieldComponents::Spectral::THREE != FieldComponents::Spectral::NOTUSED)
      {
         this->defineCoupling(FieldComponents::Spectral::THREE, CouplingInformation::TRIVIAL, 0, true, false, false, false);
         this->setExplicitTiming(FieldComponents::Spectral::THREE, ExplicitTiming::LINEAR);
      }
   }

   void ShellExactVectorState::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id compId) const
   {
      ShellExactStateIds::Id typeId = this->mTypeId;

      // Initialize to zero
      rNLComp.rData().setConstant(0);

      if(typeId == ShellExactStateIds::CONSTANT)
      {
         rNLComp.rData().setConstant(42);
      } else if(typeId == ShellExactStateIds::HARMONIC)
      {
         throw Exception("HARMONIC state is not implemented for vector states");

      } else if(typeId == ShellExactStateIds::TOROIDAL || typeId == ShellExactStateIds::POLOIDAL || typeId == ShellExactStateIds::TORPOL)
      {
         int nR = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
         int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
         int nPh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

         MHDFloat ro = this->eqParams().nd(NonDimensional::RO);
         MHDFloat rratio = this->eqParams().nd(NonDimensional::RRATIO);
         Array rGrid = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nR, ro, rratio);
         Array thGrid = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
         Array phGrid = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nPh);

         typedef std::vector<HarmonicModeType>::const_iterator ModeIt;
         std::pair<ModeIt, ModeIt>  modeRange;

         if(typeId == ShellExactStateIds::TOROIDAL || typeId == ShellExactStateIds::TORPOL)
         {
            Array funcSH(nPh);
            MHDFloat funcR = 1.0;

            modeRange.first = this->mSHModes.find(FieldComponents::Spectral::TOR)->second.begin();
            modeRange.second = this->mSHModes.find(FieldComponents::Spectral::TOR)->second.end();

            MHDFloat r;
            MHDFloat theta;
            nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
            for(int iR = 0; iR < nR; ++iR)
            {
               r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
               nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));

                  for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
                  {
                     int l = std::tr1::get<0>(*it);
                     int m = std::tr1::get<1>(*it);
                     MHDFloat amplitude = std::abs(std::tr1::get<2>(*it));

                     if(l == 0 and m == 0)
                     {
                        this->computeTor00(funcSH, compId, r, theta, phGrid);
                     } else if(l == 1 and m == 0)
                     {
                        this->computeTor10(funcSH, compId, r, theta, phGrid);
                     } else if(l == 1 and m == 1)
                     {
                        this->computeTor11(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 0)
                     {
                        this->computeTor20(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 1)
                     {
                        this->computeTor21(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 2)
                     {
                        this->computeTor22(funcSH, compId, r, theta, phGrid);
                     } else if(l == 5 and m == 4)
                     {
                        this->computeTor54(funcSH, compId, r, theta, phGrid);
                     } else
                     {
                        throw Exception("Requested unimplemented toroidal harmonic");
                     }
                     
                     rNLComp.addProfile(amplitude*funcR*funcSH,iTh,iR);
                  }
               }
            }
         }

         if(typeId == ShellExactStateIds::POLOIDAL || typeId == ShellExactStateIds::TORPOL)
         {
            Array funcSH(nPh);
            MHDFloat funcR = 1.0;

            modeRange.first = this->mSHModes.find(FieldComponents::Spectral::POL)->second.begin();
            modeRange.second = this->mSHModes.find(FieldComponents::Spectral::POL)->second.end();

            MHDFloat r;
            MHDFloat theta;
            nR = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
            for(int iR = 0; iR < nR; ++iR)
            {
               r = rGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iR));
               nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iR);
               for(int iTh = 0; iTh < nTh; ++iTh)
               {
                  theta = thGrid(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh, iR));
                  for(ModeIt it = modeRange.first; it != modeRange.second; ++it)
                  {
                     int l = std::tr1::get<0>(*it);
                     int m = std::tr1::get<1>(*it);
                     MHDFloat amplitude = std::abs(std::tr1::get<2>(*it));

                     if(l == 0 and m == 0)
                     {
                        this->computePol00(funcSH, compId, r, theta, phGrid);
                     } else if(l == 1 and m == 0)
                     {
                        this->computePol10(funcSH, compId, r, theta, phGrid);
                     } else if(l == 1 and m == 1)
                     {
                        this->computePol11(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 0)
                     {
                        this->computePol20(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 1)
                     {
                        this->computePol21(funcSH, compId, r, theta, phGrid);
                     } else if(l == 2 and m == 2)
                     {
                        this->computePol22(funcSH, compId, r, theta, phGrid);
                     } else if(l == 4 and m == 3)
                     {
                        this->computePol43(funcSH, compId, r, theta, phGrid);
                     } else
                     {
                        throw Exception("Requested unimplemented toroidal harmonic");
                     }

                     rNLComp.addProfile(amplitude*funcR*funcSH,iTh,iR);
                  }
               }
            }
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

    Datatypes::SpectralScalarType::PointType ShellExactVectorState::sourceTerm(FieldComponents::Spectral::Id compId, const int i, const int j, const int k) const
    {
      return Datatypes::SpectralScalarType::PointType(0);
    }

   void ShellExactVectorState::setRequirements()
   {
      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add unknown to requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(false, true, false, false));
   }

   void ShellExactVectorState::computeTor00(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor10(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= std::sqrt(3.0/Math::PI)*std::sin(theta);
         rField = Array::Ones(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor11(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= std::sqrt(3.0/(2.0*Math::PI));
         rField = phi.array().cos() + phi.array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= std::sqrt(3.0/(2.0*Math::PI))*std::cos(theta);
         rField = phi.array().cos() - phi.array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor20(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= 3.0*std::sqrt(5.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = Array::Ones(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor21(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta);
         rField = phi.array().cos() + phi.array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= std::sqrt(15.0/(2.0*Math::PI))*std::cos(2.0*theta);
         rField = phi.array().cos() - phi.array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor22(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= -std::sqrt(15.0/(2.0*Math::PI))*std::sin(theta);
         rField = (2.0*phi).array().cos() + (2.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= -std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta)*std::sin(theta);
         rField = (2.0*phi).array().cos() - (2.0*phi).array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computeTor54(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Toroidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= -(3.0/2.0)*std::sqrt(385.0/(2.0*Math::PI))*std::cos(theta)*std::pow(std::sin(theta),3);
         rField = (4.0*phi).array().cos() + (4.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0;
         factor *= -(3.0/16.0)*std::sqrt(385.0/(2.0*Math::PI))*(3.0 + 5.0*std::cos(2.0*theta))*std::pow(std::sin(theta),3);
         rField = (4.0*phi).array().cos() - (4.0*phi).array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol00(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol10(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= 2.0*std::sqrt(3.0/Math::PI)*std::cos(theta);
         rField = Array::Ones(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -std::sqrt(3.0/Math::PI)*std::sin(theta);
         rField = Array::Ones(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol11(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= -std::sqrt(6.0/Math::PI)*std::sin(theta);
         rField = phi.array().cos() - phi.array().sin();

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -std::sqrt(3.0/(2.0*Math::PI))*std::cos(theta);
         rField = phi.array().cos() +-phi.array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= std::sqrt(3.0/(2.0*Math::PI));
         rField = phi.array().cos() + phi.array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol20(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= (3.0/2.0)*std::sqrt(5.0/Math::PI)*(1.0 + 3.0*std::cos(2.0*theta));
         rField = Array::Ones(phi.size());

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -3.0*std::sqrt(5.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = Array::Ones(phi.size());

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = 0.0;
         rField = Array::Zero(phi.size());
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol21(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= -3.0*std::sqrt(30.0/Math::PI)*std::cos(theta)*std::sin(theta);
         rField = phi.array().cos() - phi.array().sin();

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -std::sqrt(15.0/(2.0*Math::PI))*std::cos(2.0*theta);
         rField = phi.array().cos() - phi.array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta);
         rField = phi.array().cos() + phi.array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol22(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= 3.0*std::sqrt(15.0/(2.0*Math::PI))*std::pow(std::sin(theta),2);
         rField = (2.0*phi).array().cos() - (2.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= std::sqrt(15.0/(2.0*Math::PI))*std::cos(theta)*std::sin(theta);
         rField = (2.0*phi).array().cos() - (2.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -std::sqrt(15.0/(2.0*Math::PI))*std::sin(theta);
         rField = (2.0*phi).array().cos() + (2.0*phi).array().sin();
      }
      rField *= factor;
   }

   void ShellExactVectorState::computePol43(Array& rField, FieldComponents::Physical::Id compId, const MHDFloat r, const MHDFloat theta, const Array& phi) const
   {
      // Poloidal part
      MHDFloat factor = 1.0;
      if(compId == FieldComponents::Physical::R)
      {
         factor = (std::pow(r,4) + std::pow(r,3) + std::pow(r,2) + r + 1.0)/r;
         factor *= -15.0*std::sqrt(35.0/Math::PI)*std::cos(theta)*std::pow(std::sin(theta),3);
         rField = (3.0*phi).array().cos() - (3.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::THETA)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= -(3.0/4.0)*std::sqrt(35.0/Math::PI)*(1.0 + 2.0*std::cos(2.0*theta))*std::pow(std::sin(theta),2);
         rField = (3.0*phi).array().cos() - (3.0*phi).array().sin();

      } else if(compId == FieldComponents::Physical::PHI)
      {
         factor = (5.0*std::pow(r,4) + 4.0*std::pow(r,3) + 3.0*std::pow(r,2) + 2.0*r + 1.0)/r;
         factor *= (9.0/4.0)*std::sqrt(35.0/Math::PI)*std::cos(theta)*std::pow(std::sin(theta),2);
         rField = (3.0*phi).array().cos() + (3.0*phi).array().sin();
      }
      rField *= factor;
   }

}
}
