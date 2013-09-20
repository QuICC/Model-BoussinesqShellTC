/** 
 * @file TestSpatialSchemeForwardScalar.cpp
 * @brief Source of the implementation of a test equation for the TTT scheme with exact known scalar physical space solution
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
#include "Equations/Tests/TestSpatialSchemeForwardScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestSpatialSchemeForwardScalar::TestSpatialSchemeForwardScalar()
      : IScalarEquation(SharedEquationParameters(new EquationParameters())), mTypeId(CONSTANT)
   {
   }

   TestSpatialSchemeForwardScalar::~TestSpatialSchemeForwardScalar()
   {
   }

   void TestSpatialSchemeForwardScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestSpatialSchemeForwardScalar::setSolutionType(const TestSpatialSchemeForwardScalar::SolutionTypeId id)
   {
      this->mTypeId = id;
   }

   void TestSpatialSchemeForwardScalar::setCoupling()
   {
      // Get 1D dimension (fast)
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get 2D dimension (slow)
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);
      // Get 3D dimension (medium)
      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();

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
      ArrayI blockNs(nK);
      blockNs.setConstant(nI*nJ);
      ArrayI rhsCols(nK);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nK, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void TestSpatialSchemeForwardScalar::setRequirements()
   {
      // Set simple requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, true));
   }

   void TestSpatialSchemeForwardScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat k_;
      MHDFloat j_;
      MHDFloat i_;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK));
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);
               rNLComp.setPoint(this->scalarPoint(i_, j_, k_), iI, iJ, iK);
            }
         }
      }
   }

   MHDFloat TestSpatialSchemeForwardScalar::computeScalarError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat k_;
      MHDFloat j_;
      MHDFloat i_;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK));
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);
               error = std::max(error, std::abs(this->unknown().dom(0).phys().point(iI,iJ,iK) - this->scalarPoint(i_, j_, k_)));
            }
         }
      }

      return error;
   }

   MHDFloat TestSpatialSchemeForwardScalar::computeGradientError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nK = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::PHYSICAL);
      int nJ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array gK = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nK);
      Array gJ = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nJ);
      Array gI = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nI);

      MHDFloat k_;
      MHDFloat j_;
      MHDFloat i_;
      nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      for(int iK = 0; iK < nK; ++iK)
      {
         k_ = gK(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iK));
         nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            j_ = gJ(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iJ, iK));
            for(int iI = 0; iI < nI; ++iI)
            {
               i_ = gI(iI);
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::ONE).point(iI,iJ,iK) - this->gradientPoint(i_, j_, k_, FieldComponents::Physical::ONE)));
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::TWO).point(iI,iJ,iK) - this->gradientPoint(i_, j_, k_, FieldComponents::Physical::TWO)));
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::THREE).point(iI,iJ,iK) - this->gradientPoint(i_, j_, k_, FieldComponents::Physical::THREE)));
            }
         }
      }

      return error;
   }

// Set test problem for TTT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TTT
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat z, const MHDFloat y, const MHDFloat x) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
         val *= (1.0 - 0.4*y + 1.0*std::pow(y,2) + 2.0*std::pow(y,4));
         val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat z, const MHDFloat y, const MHDFloat x, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
            val *= (1.0 - 0.4*y + 1.0*std::pow(y,2) + 2.0*std::pow(y,4));
            val *= (0.3 + 9.0*std::pow(x,2) - 2.8*std::pow(x,3));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
            val *= (0.4 + 2.0*y + 8.0*std::pow(y,3));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (1.0 + 4.0*z - 6.0*z);
            val *= (1.0 - 0.4*y + 1.0*std::pow(y,2) + 2.0*std::pow(y,4));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Set test problem for TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
         val *= (std::cos(5.0*th) + std::sin(7.0*th));
         val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (0.3 + 9.0*std::pow(x,2) - 2.8*std::pow(x,3));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,2));
            val *= (-5.0*std::sin(5.0*th) + 7.0*std::cos(7.0*th));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (1.0 + 4.0*z - 6.0*z);
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Set test problem for TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat x) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
         val *= (std::cos(5.0*th) + std::sin(7.0*th));
         val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat x, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (0.3 + 9.0*std::pow(x,2) - 2.8*std::pow(x,3));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
            val *= (-5.0*std::sin(5.0*th) + 7.0*std::cos(7.0*th));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (-0.9*std::sin(3.0*ph) + 10.0*std::cos(5.0*ph) + -2.8*std::sin(4.0*ph) - 6.0*std::cos(6.0*ph) + 24.0*std::sin(8.0*ph) + 0.9*std::cos(9.0*ph));
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (0.3*x + 3.0*std::pow(x,3) - 0.7*std::pow(x,4));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Set test problem for FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat kh) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
         val *= (std::cos(5.0*th) + std::sin(7.0*th));
         val *= (std::cos(3.0*kh) + std::sin(4.0*kh) + std::cos(2.0*kh));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat kh, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (-3.0*std::sin(3.0*kh) + 4.0*std::cos(4.0*kh) - 2.0*std::sin(2.0*kh));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (0.3*std::cos(3.0*ph) + 2.0*std::sin(5.0*ph) + 0.7*std::cos(4.0*ph) - 1.0*std::sin(6.0*ph) - 3.0*std::cos(8.0*ph) + 0.1*std::sin(9.0*ph));
            val *= (-5.0*std::sin(5.0*th) + 7.0*std::cos(7.0*th));
            val *= (std::cos(3.0*kh) + std::sin(4.0*kh) + std::cos(2.0*kh));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (-0.9*std::sin(3.0*ph) + 10.0*std::cos(5.0*ph) - 2.8*std::sin(4.0*ph) - 6.0*std::cos(6.0*ph) + 24.0*std::sin(8.0*ph) + 0.9*std::cos(9.0*ph));
            val *= (std::cos(5.0*th) + std::sin(7.0*th));
            val *= (std::cos(3.0*kh) + std::sin(4.0*kh) + std::cos(2.0*kh));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Set test problem for CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat r) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
         val *= (std::sin(5.0*th) + std::cos(7.0*th));
         val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat r, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4 + 6.0*r - 3.5*std::pow(r,4));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
            val *= (5.0*std::sin(5.0*th) - 7.0*std::sin(7.0*th));
            val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (1.0 + 4.0*z - 9.0*std::pow(z,2));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Set test problem for SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat r) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
         val *= (std::sin(5.0*th) + std::cos(7.0*th));
         val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat r, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4 + 6.0*r - 3.5*std::pow(r,4) - 1.2*std::pow(r,5) + 27.9*std::pow(r,8));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
            val *= (5.0*std::cos(5.0*th) - 7.0*std::sin(7.0*th));
            val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (std::cos(ph) + 3.0*std::sin(3.0*ph) + 4.0*std::cos(4.0*ph) - 5.0*std::sin(5.0*ph) - 5.0*std::cos(5.0*ph) - 7.0*std::sin(7.0*ph));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Set test problem for WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat r) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
         val *= (std::sin(5.0*th) + std::cos(7.0*th));
         val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat r, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4 + 6.0*r - 3.5*std::pow(r,4));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (-5.0 + 1.0*z + 2.0*std::pow(z,2) - 3.0*std::pow(z,3));
            val *= (5.0*std::cos(5.0*th) - 7.0*std::sin(7.0*th));
            val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (1.0 + 4.0*z - 9.0*std::pow(z,2));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Set test problem for WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF
   MHDFloat TestSpatialSchemeForwardScalar::scalarPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat r) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         MHDFloat val;
         val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
         val *= (std::sin(5.0*th) + std::cos(7.0*th));
         val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
         return val;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }

   MHDFloat TestSpatialSchemeForwardScalar::gradientPoint(const MHDFloat ph, const MHDFloat th, const MHDFloat r, const FieldComponents::Physical::Id compId) const
   {
      if(this->mTypeId == ZERO)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == CONSTANT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            MHDFloat val;
            val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (0.4 + 6.0*r - 3.5*std::pow(r,4) - 1.2*std::pow(r,5) + 27.9*std::pow(r,8));
            return val;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            MHDFloat val;
            val = (std::sin(ph) - std::cos(3.0*ph) + std::sin(4.0*ph) + std::cos(5.0*ph) - std::sin(5.0*ph) + std::cos(7.0*ph));
            val *= (5.0*std::cos(5.0*th) - 7.0*std::sin(7.0*th));
            val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
            return val;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            MHDFloat val;
            val = (std::cos(ph) + 3.0*std::sin(3.0*ph) + 4.0*std::cos(4.0*ph) - 5.0*std::sin(5.0*ph) - 5.0*std::cos(5.0*ph) - 7.0*std::sin(7.0*ph));
            val *= (std::sin(5.0*th) + std::cos(7.0*th));
            val *= (-0.4 + 0.4*r + 3.0*std::pow(r,2) - 0.7*std::pow(r,5) - 0.2*std::pow(r,6) + 3.1*std::pow(r,9));
            return val;
         }
      }
      // Unknown setup
      throw Exception("Unknown exact state");

   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

}
}
