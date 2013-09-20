/** 
 * @file TestTFTForwardScalar.cpp
 * @brief Source of the implementation of a test equation for the TFT scheme with exact known scalar physical space solution
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
#include "Equations/Tests/TestTFTForwardScalar.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   TestTFTForwardScalar::TestTFTForwardScalar()
      : IScalarEquation(SharedEquationParameters(new EquationParameters())), mTypeId(CONSTANT)
   {
   }

   TestTFTForwardScalar::~TestTFTForwardScalar()
   {
   }

   void TestTFTForwardScalar::setIdentity(const PhysicalNames::Id name)
   {
      // Set the name
      this->setName(name);

      // Set the variable requirements
      this->setRequirements();
   }

   void TestTFTForwardScalar::setSolutionType(const TestTFTForwardScalar::SolutionTypeId id)
   {
      this->mTypeId = id;
   }

   void TestTFTForwardScalar::setCoupling()
   {
      // Get X dimension
      int nX = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL);
      // Get Y dimension
      int nY = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>();
      // Get Z dimension
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::SPECTRAL);

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
      ArrayI blockNs(nY);
      blockNs.setConstant(nX*nZ);
      ArrayI rhsCols(nY);
      rhsCols.setConstant(1);
      infoIt.first->second.setSizes(nY, blockNs, rhsCols); 

      // Sort implicit fields
      infoIt.first->second.sortImplicitFields(eqId.first, FieldComponents::Spectral::SCALAR);
   }

   void TestTFTForwardScalar::setRequirements()
   {
      // Set simple requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(this->name(), FieldRequirement(true, true, true, true));
   }

   MHDFloat TestTFTForwardScalar::scalarPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x) const
   {
      if(this->mTypeId == ZERO)
      {
         return 0.0;
      } else if(this->mTypeId == CONSTANT)
      {
         return 42.0;
      } else if(this->mTypeId == EXACT)
      {
         return std::pow(z,11);
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

   MHDFloat TestTFTForwardScalar::gradientPoint(const MHDFloat z, const MHDFloat th, const MHDFloat x, const FieldComponents::Physical::Id compId) const
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
         } else
         {
            throw Exception("Unknown field component");
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
         } else
         {
            throw Exception("Unknown field component");
         }
      } else if(this->mTypeId == EXACT)
      {
         if(compId == FieldComponents::Physical::ONE)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::TWO)
         {
            return 0.0;
         } else if(compId == FieldComponents::Physical::THREE)
         {
            return 11.0*std::pow(z,10);
         } else
         {
            throw Exception("Unknown field component");
         }
      } else
      {
         throw Exception("Unknown exact state");
      }
   }

   void TestTFTForwardScalar::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array x = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
      Array th = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      Array z = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

      MHDFloat x_;
      MHDFloat th_;
      MHDFloat z_;
      for(int iX = 0; iX < nX; ++iX)
      {
         x_ = x(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iX));
         nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iX);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            th_ = th(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh));
            for(int iZ = 0; iZ < nZ; ++iZ)
            {
               z_ = z(iZ);
               rNLComp.setPoint(this->scalarPoint(z_, th_, x_), iZ, iTh, iX);
            }
         }
      }
   }

   MHDFloat TestTFTForwardScalar::computeScalarError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array x = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
      Array th = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      Array z = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

      MHDFloat x_;
      MHDFloat th_;
      MHDFloat z_;
      for(int iX = 0; iX < nX; ++iX)
      {
         x_ = x(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iX));
         nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iX);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            th_ = th(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh));
            for(int iZ = 0; iZ < nZ; ++iZ)
            {
               z_ = z(iZ);
               error = std::max(error, std::abs(this->unknown().dom(0).phys().point(iZ,iTh,iX) - this->scalarPoint(z_, th_, x_)));
            }
         }
      }

      return error;
   }

   MHDFloat TestTFTForwardScalar::computeGradientError() const
   {
      // Maximum error
      MHDFloat error = 0.0;

      int nX = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>();
      int nTh = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::PHYSICAL);
      int nZ = this->unknown().dom(0).spRes()->sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::PHYSICAL);

      Array x = Transform::TransformSelector<Dimensions::Transform::TRA1D>::Type::generateGrid(nX);
      Array th = Transform::TransformSelector<Dimensions::Transform::TRA2D>::Type::generateGrid(nTh);
      Array z = Transform::TransformSelector<Dimensions::Transform::TRA3D>::Type::generateGrid(nZ);

      MHDFloat x_;
      MHDFloat th_;
      MHDFloat z_;
      for(int iX = 0; iX < nX; ++iX)
      {
         x_ = x(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(iX));
         nTh = this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT2D>(iX);
         for(int iTh = 0; iTh < nTh; ++iTh)
         {
            th_ = th(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT2D>(iTh));
            for(int iZ = 0; iZ < nZ; ++iZ)
            {
               z_ = z(iZ);
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::ONE).point(iZ,iTh,iX) - this->gradientPoint(z_, th_, x_, FieldComponents::Physical::ONE)));
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::TWO).point(iZ,iTh,iX) - this->gradientPoint(z_, th_, x_, FieldComponents::Physical::TWO)));
               error = std::max(error, std::abs(this->unknown().dom(0).grad().comp(FieldComponents::Physical::THREE).point(iZ,iTh,iX) - this->gradientPoint(z_, th_, x_, FieldComponents::Physical::THREE)));
            }
         }
      }

      return error;
   }

}
}
