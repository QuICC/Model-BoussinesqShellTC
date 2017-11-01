/** 
 * @file SpatialSchemeBackwardScalar.cpp
 * @brief Source of the implementation of a test equation for the spatial schemes with exact known scalar spectral space solution
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

/// Define small macros allowing to convert to string
#define MAKE_STR_X( _P ) # _P
#define MAKE_STR( _P ) MAKE_STR_X( _P )

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Test/SpatialSchemeBackwardScalar.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace QuICC {

namespace Equations {

namespace Test {

   SpatialSchemeBackwardScalar::SpatialSchemeBackwardScalar()
      : IScalarEquation(SharedEquationParameters(new EquationParameters())), mTypeId(CONSTANT)
   {
   }

   SpatialSchemeBackwardScalar::~SpatialSchemeBackwardScalar()
   {
   }

   void SpatialSchemeBackwardScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void SpatialSchemeBackwardScalar::setSolutionType(const SpatialSchemeBackwardScalar::SolutionTypeId id)
   {
      this->mTypeId = id;
   }

   void SpatialSchemeBackwardScalar::setCoupling()
   {
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void SpatialSchemeBackwardScalar::setRequirements()
   {
      // Set simple requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, false));
   }

   void SpatialSchemeBackwardScalar::init()
   {
      // Call parent implentation
      IEquation::init();
      
      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM);

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

   void SpatialSchemeBackwardScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
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

   MHDFloat SpatialSchemeBackwardScalar::computeScalarError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nK = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      int nI = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D,Dimensions::Space::TRANSFORM);

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
#ifdef QUICC_SPATIALSCHEME_TTT
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = 2.0*(k+1)*std::pow(-1,i);
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val = 1.0;

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      } 

      return val;
   }
#endif //QUICC_SPATIALSCHEME_TTT

// Set test problem for TFT scheme
#ifdef QUICC_SPATIALSCHEME_TFT
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
          val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            if(k == 0)
            {
               val = MHDComplex(2.0*(k+1)*std::pow(-1,i), 0.0);
            } else
            {
               val = MHDComplex(2.0*(k+1)*std::pow(-1,i), 1.5*(k+1)*std::pow(-1,i));
            }
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         if(k == 0)
         {
            val = MHDComplex(1.0,0.0);
         } else
         {
            val = MHDComplex(1.0,1.0);
         }

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_TFT

// Set test problem for TFF scheme
#ifdef QUICC_SPATIALSCHEME_TFF
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val =  0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            if(j > 0)
            {
               val = MHDComplex(2.0*(k+1)*std::pow(-1,i), 1.5*(k+1)*std::pow(-1,i));
            } else
            {
              val = 1.0;
            }
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         if(j > 0)
         {
            val = MHDComplex(1.0,1.0);
         } else
         {
            val = 1.0;
         }

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_TFF

// Set test problem for FFF scheme
#ifdef QUICC_SPATIALSCHEME_FFF
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = MHDCompleX(2.0*(k+1)*std::pow(-1,i),1.5*(k+1)*std::pow(-1,i));
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val = MHDComplex(1.0,1.0);
      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_FFF

// Set test problem for CFT scheme
#ifdef QUICC_SPATIALSCHEME_CFT
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = MHDComplex(2.0*(k+1)*std::pow(-1,i),1.5*(k+1)*std::pow(-1,i));
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val = MHDComplex(1.0,1.0);

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      } 

      return val;
   }
#endif //QUICC_SPATIALSCHEME_CFT

// Set test problem for SLF scheme
#ifdef QUICC_SPATIALSCHEME_SLF
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val 42.0;
         } else
         {
            val 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = MHDComplex(2.0*(k+1)*std::pow(-1,i),1.5*(k+1)*std::pow(-1,i));
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val =  MHDComplex(1.0,1.0);

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_SLF

// Set test problem for WFT scheme
#ifdef QUICC_SPATIALSCHEME_WFT
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = MHDComplex(2.0*(k+1)*std::pow(-1,i),1.5*(k+1)*std::pow(-1,i));
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val = MHDComplex(1.0,1.0);

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_WFT

// Set test problem for WLF scheme
#ifdef QUICC_SPATIALSCHEME_WLF
   Datatypes::SpectralScalarType::PointType SpatialSchemeBackwardScalar::scalarPoint(const int i, const int j, const int k) const
   {
      Datatypes::SpectralScalarType::PointType val;

      if(this->mTypeId == ZERO)
      {
         val = 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         if(i == 0 && j == 0 && k == 0)
         {
            val = 42.0;
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == EXACT)
      {
         if(i < 10 && j < 10 && k < 10)
         {
            val = MHDComplex(2.0*(k+1)*std::pow(-1,i),1.5*(k+1)*std::pow(-1,i));
         } else
         {
            val = 0.0;
         }
      } else if(this->mTypeId == FULL)
      {
         val = MHDComplex(1.0,1.0);

      // Unknown setup
      } else
      {
         throw Exception("Unknown exact state");
      }

      return val;
   }
#endif //QUICC_SPATIALSCHEME_WLF

}
}
}