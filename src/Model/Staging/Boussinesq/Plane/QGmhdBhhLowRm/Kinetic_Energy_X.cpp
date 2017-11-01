/** 
 * @file Kinetic_Energy_X.cpp
 * @brief Source of the implementation of the local thermal dissipation computation in the F-plane 3DQG model
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

// Class include
//
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/QGmhdBhhLowRm/Kinetic_Energy_X.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/FieldIdsTools.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace QGmhdBhhLowRm {

   Kinetic_Energy_X::Kinetic_Energy_X(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   Kinetic_Energy_X::~Kinetic_Energy_X()
   {
   }

   void Kinetic_Energy_X::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 1, true, true);
   }

   void Kinetic_Energy_X::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);


      /// 
      /// Computation:
      ///   \f$ 0.5*VELOCITYX^2  \f$
      ///
      
      rNLComp.setData(0.5*(this->scalar(PhysicalNames::VELOCITYX).dom(0).phys().data().array()*this->scalar(PhysicalNames::VELOCITYX).dom(0).phys().data().array()).matrix());
}

   Datatypes::SpectralScalarType::PointType Kinetic_Energy_X::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void Kinetic_Energy_X::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::KINETIC_ENERGY_X);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add Kinetic_Energy_X requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::KINETIC_ENERGY_X, FieldRequirement(true, true, true, false, false, false));

      // Add VELOCITYX requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::VELOCITYX, FieldRequirement(true, true, true, false, false, false));


   }
}
}
}
}
}
