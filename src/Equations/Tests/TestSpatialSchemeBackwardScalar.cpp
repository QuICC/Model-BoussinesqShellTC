/** 
 * @file TestSpatialSchemeBackwardScalar.cpp
 * @brief Source of the implementation of a test equation for the spatial schemes with exact known scalar spectral space solution
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
#include "Equations/Tests/TestSpatialSchemeBackwardScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

#include <iostream>
namespace GeoMHDiSCC {

namespace Equations {

   TestSpatialSchemeBackwardScalar::TestSpatialSchemeBackwardScalar()
      : IScalarEquation(SharedEquationParameters(new EquationParameters())), mTypeId(CONSTANT)
   {
   }

   TestSpatialSchemeBackwardScalar::~TestSpatialSchemeBackwardScalar()
   {
   }

   void TestSpatialSchemeBackwardScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestSpatialSchemeBackwardScalar::setSolutionType(const TestSpatialSchemeBackwardScalar::SolutionTypeId id)
   {
      this->mTypeId = id;
   }

   void TestSpatialSchemeBackwardScalar::setCoupling()
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

   void TestSpatialSchemeBackwardScalar::setRequirements()
   {
      // Set simple requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   void TestSpatialSchemeBackwardScalar::init()
   {
      // Call parent implentation
      IEquation::init();
      
      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL);

      int iI_;
      int iJ_;
      int iK_;
      for(int iK = 0; iK < nK; ++iK)
      {
         iK_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK);
         int nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            iJ_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ, iK);
            for(int iI = 0; iI < nI; ++iI)
            {
               iI_ = iI;
               this->rUnknown().rDom(0).rPerturbation().setPoint(this->scalarPoint(iI_, iJ_, iK_),iI,iJ,iK);
            }
         }
      }
   }

   void TestSpatialSchemeBackwardScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      for(int iK = 0; iK < nK; ++iK)
      {
         int nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            for(int iI = 0; iI < nI; ++iI)
            {
               rNLComp.setProfile(this->unknown().dom(0).phys().profile(iJ,iK), iJ, iK);
            }
         }
      }
   }

   MHDFloat TestSpatialSchemeBackwardScalar::computeScalarError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::SPECTRAL);

      int iI_;
      int iJ_;
      int iK_;
      for(int iK = 0; iK < nK; ++iK)
      {
         iK_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iK);
         int nJ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(iK);
         for(int iJ = 0; iJ < nJ; ++iJ)
         {
            iJ_ = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iJ, iK);
            for(int iI = 0; iI < nI; ++iI)
            {
               iI_ = iI;
               error = std::max(error, std::abs(this->unknown().dom(0).perturbation().point(iI,iJ,iK) - this->scalarPoint(iI_, iJ_, iK_)));
            }
         }
      }

      return error;
   }

// Set test problem for TTT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TTT
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TTT

// Set test problem for TFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFT
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFT

// Set test problem for TFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_TFF
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_TFF

// Set test problem for FFF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_FFF
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_FFF

// Set test problem for CFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_CFT
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_CFT

// Set test problem for SLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_SLF
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_SLF

// Set test problem for WFT scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WFT
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WFT

// Set test problem for WLF scheme
#ifdef GEOMHDISCC_SPATIALSCHEME_WLF
   MHDFloat TestSpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            return 42.0;
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            return 2.0*k*std::pow(-1,i);
         } else
         {
            return 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         return 1.0;
      }

      // Unknown setup
      throw Exception("Unknown exact state");
   }
#endif //GEOMHDISCC_SPATIALSCHEME_WLF

}
}
