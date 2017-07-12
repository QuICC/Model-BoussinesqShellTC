/** 
 * @file HumanToId.cpp
 * @brief Source of human strings to enum id converters
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */
// Configuration includes
//
#include "Exceptions/Exception.hpp"

// Configuration includes
//

// System includes
//

// External includes
//

// Class include
//
#include "IoTools/HumanToId.hpp"

// Project includes
//
#include "IoTools/IdToHuman.hpp"

namespace QuICC {

namespace IoTools {

   NonDimensional::Id HumanToId::toNd(const std::string& id)
   {
      //
      // Nondimensional numbers
      // 
      if(id == IdToHuman::toTag(NonDimensional::EADY))
      {
         return NonDimensional::EADY;

      } else if(id == IdToHuman::toTag(NonDimensional::EKMAN))
      {
         return NonDimensional::EKMAN;

      } else if(id == IdToHuman::toTag(NonDimensional::MAGREYNOLDS))
      {
    	  return NonDimensional::MAGREYNOLDS;

      } else if(id == IdToHuman::toTag(NonDimensional::MAGEKMAN))
      {
         return NonDimensional::MAGEKMAN;

      } else if(id == IdToHuman::toTag(NonDimensional::MAGPRANDTL))
      {
         return NonDimensional::MAGPRANDTL;

      } else if(id == IdToHuman::toTag(NonDimensional::POINCARE))
      {
         return NonDimensional::POINCARE;

      } else if(id == IdToHuman::toTag(NonDimensional::MODELSASSER))
      {
    	  return NonDimensional::MODELSASSER;

      } else if(id == IdToHuman::toTag(NonDimensional::PRANDTL))
      {
         return NonDimensional::PRANDTL;

      } else if(id == IdToHuman::toTag(NonDimensional::RAYLEIGH))
      {
         return NonDimensional::RAYLEIGH;

      } else if(id == IdToHuman::toTag(NonDimensional::CHANDRASEKHAR))
      {
         return NonDimensional::CHANDRASEKHAR;

      } else if(id == IdToHuman::toTag(NonDimensional::ROBERTS))
      {
         return NonDimensional::ROBERTS;

      } else if(id == IdToHuman::toTag(NonDimensional::ROSSBY))
      {
         return NonDimensional::ROSSBY;

      } else if(id == IdToHuman::toTag(NonDimensional::TAYLOR))
      {
         return NonDimensional::TAYLOR;

      // Geometrical numbers
      //
      } else if(id == IdToHuman::toTag(NonDimensional::RO))
      {
         return NonDimensional::RO;

      } else if(id == IdToHuman::toTag(NonDimensional::RRATIO))
      {
         return NonDimensional::RRATIO;

      // Flags
      //
      } else if(id == IdToHuman::toTag(NonDimensional::HEATING))
      {
         return NonDimensional::HEATING;

      // Axis scaling factors
      } else if(id == IdToHuman::toTag(NonDimensional::SCALE1D))
      {
         return NonDimensional::SCALE1D;

      } else if(id == IdToHuman::toTag(NonDimensional::SCALE2D))
      {
         return NonDimensional::SCALE2D;

      } else if(id == IdToHuman::toTag(NonDimensional::SCALE3D))
      {
         return NonDimensional::SCALE3D;

      //
      // Greek alphabet
      } else if(id == IdToHuman::toTag(NonDimensional::ALPHA))
      {
         return NonDimensional::ALPHA;

      } else if(id == IdToHuman::toTag(NonDimensional::BETA))
      {
         return NonDimensional::BETA;

      } else if(id == IdToHuman::toTag(NonDimensional::GAMMA))
      {
         return NonDimensional::GAMMA;

      } else if(id == IdToHuman::toTag(NonDimensional::DELTA))
      {
         return NonDimensional::DELTA;

      } else if(id == IdToHuman::toTag(NonDimensional::EPSILON))
      {
         return NonDimensional::EPSILON;

      } else if(id == IdToHuman::toTag(NonDimensional::ZETA))
      {
         return NonDimensional::ZETA;

      } else if(id == IdToHuman::toTag(NonDimensional::ETA))
      {
         return NonDimensional::ETA;

      } else if(id == IdToHuman::toTag(NonDimensional::THETA))
      {
         return NonDimensional::THETA;

      } else if(id == IdToHuman::toTag(NonDimensional::IOTA))
      {
         return NonDimensional::IOTA;

      } else if(id == IdToHuman::toTag(NonDimensional::KAPPA))
      {
         return NonDimensional::KAPPA;

      } else if(id == IdToHuman::toTag(NonDimensional::LAMBDA))
      {
         return NonDimensional::LAMBDA;

      } else if(id == IdToHuman::toTag(NonDimensional::MU))
      {
         return NonDimensional::MU;

      } else if(id == IdToHuman::toTag(NonDimensional::NU))
      {
         return NonDimensional::NU;

      } else if(id == IdToHuman::toTag(NonDimensional::XI))
      {
         return NonDimensional::XI;

      } else if(id == IdToHuman::toTag(NonDimensional::OMICRON))
      {
         return NonDimensional::OMICRON;

      } else if(id == IdToHuman::toTag(NonDimensional::PI))
      {
         return NonDimensional::PI;

      } else if(id == IdToHuman::toTag(NonDimensional::RHO))
      {
         return NonDimensional::RHO;

      } else if(id == IdToHuman::toTag(NonDimensional::SIGMA))
      {
         return NonDimensional::SIGMA;

      } else if(id == IdToHuman::toTag(NonDimensional::TAU))
      {
         return NonDimensional::TAU;

      } else if(id == IdToHuman::toTag(NonDimensional::UPSILON))
      {
         return NonDimensional::UPSILON;

      } else if(id == IdToHuman::toTag(NonDimensional::PHI))
      {
         return NonDimensional::PHI;

      } else if(id == IdToHuman::toTag(NonDimensional::CHI))
      {
         return NonDimensional::CHI;

      } else if(id == IdToHuman::toTag(NonDimensional::PSI))
      {
         return NonDimensional::PSI;

      } else if(id == IdToHuman::toTag(NonDimensional::OMEGA))
      {
         return NonDimensional::OMEGA;

      //
      // Special flags
      } else if(id == IdToHuman::toTag(NonDimensional::ELEVATOR))
      {
         return NonDimensional::ELEVATOR;

      } else if(id == IdToHuman::toTag(NonDimensional::FAST_MEAN))
      {
         return NonDimensional::FAST_MEAN;

      } else if(id == IdToHuman::toTag(NonDimensional::RESCALED))
      {
         return NonDimensional::RESCALED;

      } else
      {
         throw Exception("Unknown string to ID conversion requested (Nondimensional)");
      }
   }

   ModelOperator::Id HumanToId::toModOp(const std::string& id)
   {
      if(id == IdToHuman::toString(ModelOperator::TIME))
      {
         return ModelOperator::TIME;

      } else if(id == IdToHuman::toString(ModelOperator::IMPLICIT_LINEAR))
      {
         return ModelOperator::IMPLICIT_LINEAR;

      } else if(id == IdToHuman::toString(ModelOperator::BOUNDARY))
      {
         return ModelOperator::BOUNDARY;

      } else if(id == IdToHuman::toString(ModelOperator::EXPLICIT_LINEAR))
      {
         return ModelOperator::EXPLICIT_LINEAR;

      } else if(id == IdToHuman::toString(ModelOperator::EXPLICIT_NONLINEAR))
      {
         return ModelOperator::EXPLICIT_NONLINEAR;

      } else if(id == IdToHuman::toString(ModelOperator::EXPLICIT_NEXTSTEP))
      {
         return ModelOperator::EXPLICIT_NEXTSTEP;

      } else if(id == IdToHuman::toString(ModelOperator::STENCIL))
      {
         return ModelOperator::STENCIL;

      } else
      {
         throw Exception("Unknown string to ID conversion requested (ModelOperator)");
      }
   }

   PhysicalNames::Id HumanToId::toPhys(const std::string& id)
   {
      if(id == IdToHuman::toTag(PhysicalNames::CODENSITY))
      {
         return PhysicalNames::CODENSITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::BX))
      {
         return PhysicalNames::BX;

      } else if(id == IdToHuman::toTag(PhysicalNames::BY))
      {
         return PhysicalNames::BY;

      } else if(id == IdToHuman::toTag(PhysicalNames::DENSITY))
      {
         return PhysicalNames::DENSITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::DZ_MEANTEMPERATURE))
      {
         return PhysicalNames::DZ_MEANTEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::DX_MEANTEMPERATURE))
      {
         return PhysicalNames::DX_MEANTEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::EMFX))
      {
         return PhysicalNames::EMFX;

      } else if(id == IdToHuman::toTag(PhysicalNames::EMFY))
      {
         return PhysicalNames::EMFY;

      } else if(id == IdToHuman::toTag(PhysicalNames::ENTROPY))
      {
         return PhysicalNames::ENTROPY;

      } else if(id == IdToHuman::toTag(PhysicalNames::FBX))
      {
         return PhysicalNames::FBX;

      } else if(id == IdToHuman::toTag(PhysicalNames::FBY))
      {
         return PhysicalNames::FBY;

      } else if(id == IdToHuman::toTag(PhysicalNames::FBZ))
      {
         return PhysicalNames::FBZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::MAGNETIC))
      {
         return PhysicalNames::MAGNETIC;

      } else if(id == IdToHuman::toTag(PhysicalNames::PRESSURE))
      {
         return PhysicalNames::PRESSURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::STREAMFUNCTION))
      {
         return PhysicalNames::STREAMFUNCTION;

      } else if(id == IdToHuman::toTag(PhysicalNames::TEMPERATURE))
      {
         return PhysicalNames::TEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::VELOCITY))
      {
         return PhysicalNames::VELOCITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::VELOCITYX))
      {
         return PhysicalNames::VELOCITYX;

      } else if(id == IdToHuman::toTag(PhysicalNames::VELOCITYY))
      {
         return PhysicalNames::VELOCITYY;

      } else if(id == IdToHuman::toTag(PhysicalNames::VELOCITYZ))
      {
         return PhysicalNames::VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::ZONAL_VELOCITY))
      {
         return PhysicalNames::ZONAL_VELOCITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::NONZONAL_VELOCITY))
      {
         return PhysicalNames::NONZONAL_VELOCITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::ROTATEDGEOSTROPHIC_VELOCITY))
	  {
    	 return PhysicalNames::ROTATEDGEOSTROPHIC_VELOCITY;

	  } else if(id == IdToHuman::toTag(PhysicalNames::VORTICITY))
	  {
         return PhysicalNames::VORTICITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::VORTICITYX))
      {
         return PhysicalNames::VORTICITYX;

      } else if(id == IdToHuman::toTag(PhysicalNames::VORTICITYY))
      {
         return PhysicalNames::VORTICITYY;

      } else if(id == IdToHuman::toTag(PhysicalNames::VORTICITYZ))
      {
         return PhysicalNames::VORTICITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::PHI))
      {
         return PhysicalNames::PHI;

      } else if(id == IdToHuman::toTag(PhysicalNames::NO_STREAMFUNCTION))
      {
         return PhysicalNames::NO_STREAMFUNCTION;

      } else if(id == IdToHuman::toTag(PhysicalNames::NO_VELOCITYZ))
      {
         return PhysicalNames::NO_VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::NO_VORTICITYZ))
      {
         return PhysicalNames::NO_VORTICITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_TEMPERATURE))
      {
         return PhysicalNames::TILTED_TEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_STREAMFUNCTION))
      {
         return PhysicalNames::TILTED_STREAMFUNCTION;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_VELOCITYZ))
      {
         return PhysicalNames::TILTED_VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_NO_STREAMFUNCTION))
      {
         return PhysicalNames::TILTED_NO_STREAMFUNCTION;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_NO_VELOCITYZ))
      {
         return PhysicalNames::TILTED_NO_VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::TILTED_NO_VORTICITYZ))
      {
         return PhysicalNames::TILTED_NO_VORTICITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_TEMPERATURE))
      {
         return PhysicalNames::FLUCT_TEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_MAGNETIC))
      {
         return PhysicalNames::FLUCT_MAGNETIC;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_MAGNETICX))
      {
         return PhysicalNames::FLUCT_MAGNETICX;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_MAGNETICY))
      {
         return PhysicalNames::FLUCT_MAGNETICY;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_MAGNETICZ))
      {
         return PhysicalNames::FLUCT_MAGNETICZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_VELOCITY))
      {
         return PhysicalNames::FLUCT_VELOCITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_VELOCITYX))
      {
         return PhysicalNames::FLUCT_VELOCITYX;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_VELOCITYY))
      {
         return PhysicalNames::FLUCT_VELOCITYY;

      } else if(id == IdToHuman::toTag(PhysicalNames::FLUCT_VELOCITYZ))
      {
         return PhysicalNames::FLUCT_VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_TEMPERATURE))
      {
         return PhysicalNames::MEAN_TEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_MAGNETIC))
      {
         return PhysicalNames::MEAN_MAGNETIC;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_MAGNETICX))
      {
         return PhysicalNames::MEAN_MAGNETICX;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_MAGNETICY))
      {
         return PhysicalNames::MEAN_MAGNETICY;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_MAGNETICZ))
      {
         return PhysicalNames::MEAN_MAGNETICZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_VELOCITY))
      {
         return PhysicalNames::MEAN_VELOCITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_VELOCITYX))
      {
         return PhysicalNames::MEAN_VELOCITYX;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_VELOCITYY))
      {
         return PhysicalNames::MEAN_VELOCITYY;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEAN_VELOCITYZ))
      {
         return PhysicalNames::MEAN_VELOCITYZ;

      } else if(id == IdToHuman::toTag(PhysicalNames::KINETIC_ENERGY))
      {
         return PhysicalNames::KINETIC_ENERGY;

      } else if(id == IdToHuman::toTag(PhysicalNames::ZONAL_KINETIC_ENERGY))
      {
         return PhysicalNames::ZONAL_KINETIC_ENERGY;

      } else if(id == IdToHuman::toTag(PhysicalNames::NONZONAL_KINETIC_ENERGY))
      {
         return PhysicalNames::NONZONAL_KINETIC_ENERGY;

      } else if(id == IdToHuman::toTag(PhysicalNames::IMPOSED_MAGNETIC))
      {
         return PhysicalNames::IMPOSED_MAGNETIC;

      } else
      {
         throw Exception("Unknown string to ID conversion requested (PhysicalNames)");
      }
   }

   FieldComponents::Spectral::Id HumanToId::toComp(const std::string& id)
   {
      if(id == IdToHuman::toTag(FieldComponents::Spectral::SCALAR))
      {
         return FieldComponents::Spectral::SCALAR;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::ONE))
      {
         return FieldComponents::Spectral::ONE;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::TWO))
      {
         return FieldComponents::Spectral::TWO;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::THREE))
      {
         return FieldComponents::Spectral::THREE;
      } else
      {
         throw Exception("Unknown string to ID conversion requested (FieldComponents::Spectral)");
      }
   }
}

}
