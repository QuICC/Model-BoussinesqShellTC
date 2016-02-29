/** 
 * @file BoussinesqDynamo3DQGfby.cpp
 * @brief Source of the implementation of the mean heat computation in the F-plane 3DQG model
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//

// System includes
//

// External includes

// Class include
//
#include "Equations/Asymptotics/FPlane3DQG/Boussinesq/BoussinesqDynamo3DQGfby.hpp"

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   BoussinesqDynamo3DQGfby::BoussinesqDynamo3DQGfby(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   BoussinesqDynamo3DQGfby::~BoussinesqDynamo3DQGfby()
   {
   }

   void BoussinesqDynamo3DQGfby::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::DIAGNOSTIC, 1, true, true);
   }

   void BoussinesqDynamo3DQGfby::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);

      // Get paramters
      MHDFloat MPr = this->eqParams().nd(NonDimensional::MAGPRANDTL);

      /// 
      /// Computation:
      ///   \f$ MPr(Bx dxy\Psi + By dyy\Psi) \f$
      ///
      rNLComp.setData((-MPr*(this->scalar(PhysicalNames::BX).dom(0).phys().data().array()*this->scalar(PhysicalNames::STREAMFUNCTION).dom(0).grad2().comp(FieldComponents::Physical::X,FieldComponents::Physical::X).data().array() + this->scalar(PhysicalNames::BY).dom(0).phys().data().array*this->scalar(PhysicalNames:STREAMFUNCTION).dom().grad2().comp(FieldComponents::Physical::X,FieldComponents:Physical::Y).data().array())).matrix());
   }

   Datatypes::SpectralScalarType::PointType BoussinesqDynamo3DQGfby::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
   {
      // Assert on scalar component is used
      assert(compId == FieldComponents::Spectral::SCALAR);

      if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(iY) == 0)
      {
         if(this->unknown().dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(iZ,iY) == 0)
         {
            if(iX == 0)
            {
               return Datatypes::SpectralScalarType::PointType(-1.0);
            }
         }
      } 

      return Datatypes::SpectralScalarType::PointType(0.0);
   }

   void BoussinesqDynamo3DQGfby::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::FBY);

      // Set solver timing
      this->setSolveTiming(SolveTiming::DIAGNOSTIC);

      // Add fBY requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::FBY, FieldRequirement(true, true, true, false));

      // Add BX requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BX, FieldRequirement(true, true, true, false));

      // Add BY requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::BY, FieldRequirement(true, true, true, false));

      // Add streamfunction requirements: is scalar?, need spectral?, need physical?, need diff? need curl? need diff2?
      this->mRequirements.addField(PhysicalNames::STREAMFUNCTION, FieldRequirement(true, true, false, false, true));
   }

//      // Restrict components of 2nd order gradient
//      // Make upper triangular matrix
      MatrixB   compsTen = MatrixB::Constant(3,3, true);
      compsTen.triangularView<Eigen::StrictlyLower>().setZero();
      // Don't compute DxDz derivative (order doesn't matter, but use fieldPairSym)
      std::pair<int,int> idx = fieldPairSym(FieldComponents::Physical::Z,FieldComponents::Physical::X);
      compsTen(idx.first, idx.second) = false;
      // Don't compute DzDy derivative (order doesn't matter, but use fieldPairSym)
      idx = fieldPairSym(FieldComponents::Physical::Y,FieldComponents::Physical::Z);
      compsTen(idx.first, idx.second) = false;
      // Don't compute DzDz derivative (order doesn't matter, but use fieldPairSym)
      std::pair<int,int> idx = fieldPairSym(FieldComponents::Physical::Z,FieldComponents::Physical::Z);
      compsTen(idx.first, idx.second) = false;
      std::map<FieldComponents::Spectral::Id,MatrixB>  grad2Comps;
      grad2Comps.insert(std::make_pair(FieldComponents::Spectral::SCALAR, compsTen));
//
      this->updateFieldRequirements(PhysicalNames::STREAMFUNCTION).updateGradient2(grad2Comps);
}
}
