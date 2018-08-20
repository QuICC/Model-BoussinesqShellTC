/** 
 * @file DissTh.cpp
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
#include MAKE_STR( QUICC_MODEL_PATH/Boussinesq/Plane/F3DQG/DissTh.hpp )

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Base/MathConstants.hpp"
#include "Enums/NonDimensional.hpp"
#include "Enums/FieldIdsTools.hpp"
#include "TypeSelectors/EquationEigenSelector.hpp"
#include <cmath>

namespace QuICC {

namespace Equations {

namespace Boussinesq {

namespace Plane {

namespace F3DQG {

   DissTh::DissTh(SharedEquationParameters spEqParams)
      : IScalarEquation(spEqParams)
   {
      // Set the variable requirements
      this->setRequirements();
   }

   DissTh::~DissTh()
   {
   }

   void DissTh::setCoupling()
   {	
      // 1: want index to start at 1 because of inverse laplacian, T, T?
      this->defineCoupling(FieldComponents::Spectral::SCALAR, CouplingInformation::TRIVIAL, 0, true, false);
   }

   void DissTh::computeNonlinear(Datatypes::PhysicalScalarType& rNLComp, FieldComponents::Physical::Id id) const
   {
      // Assert on scalar component is used
      assert(id == FieldComponents::Physical::SCALAR);


      /// 
      /// Computation:
      ///   \f$ DZ_MEANTEMPERATURE^2 + (\nabla TEMPERATURE)^2  \f$
      ///
      /// A +1 is added to solve for the -1 level in the background always present in the field. Not sure where this comes from, but probably is a m=0 problem. Or a scaling transformation gone wrong. 
	rNLComp.setData((1+pow(this->scalar(PhysicalNames::DZ_MEANTEMPERATURE).dom(0).phys().data().array(),2) + pow(this->scalar(PhysicalNames::TEMPERATURE).dom(0).grad().comp(FieldComponents::Physical::X).data().array(),2) + pow(this->scalar(PhysicalNames::TEMPERATURE).dom(0).grad().comp(FieldComponents::Physical::Y).data().array(),2)).matrix());

}

   Datatypes::SpectralScalarType::PointType DissTh::sourceTerm(FieldComponents::Spectral::Id compId, const int iX, const int iZ, const int iY) const
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

   void DissTh::setRequirements()
   {
      // Set fluctuating bx field as equation unknown
      this->setName(PhysicalNames::DISSTH);

      // Set solver timing
      this->setSolveTiming(SolveTiming::AFTER);

      // Add DissTh requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DISSTH, FieldRequirement(true, true, true, false, false, false));

      // Add DZ_MEANTEMPERATURE requirements: is scalar?, need spectral?, need physical?, need diff?
      this->mRequirements.addField(PhysicalNames::DZ_MEANTEMPERATURE, FieldRequirement(true, true, true, false, false, false));


      // Add TEMPERATURE requirements: is scalar?, need spectral?, need physical?, need diff? need curl? need diff2?
      this->mRequirements.addField(PhysicalNames::TEMPERATURE, FieldRequirement(true, true, true, true, false, false));

   }
}
}
}
}
}
