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

namespace GeoMHDiSCC {

namespace IoTools {

   NonDimensional::Id HumanToId::toNd(const std::string& id)
   {
      if(id == IdToHuman::toTag(NonDimensional::CHI))
      {
         return NonDimensional::CHI;

      } else if(id == IdToHuman::toTag(NonDimensional::EKMAN))
      {
         return NonDimensional::EKMAN;

      } else if(id == IdToHuman::toTag(NonDimensional::GAMMA))
      {
         return NonDimensional::GAMMA;

      } else if(id == IdToHuman::toTag(NonDimensional::MAGEKMAN))
      {
         return NonDimensional::MAGEKMAN;

      } else if(id == IdToHuman::toTag(NonDimensional::MAGPRANDTL))
      {
         return NonDimensional::MAGPRANDTL;

      } else if(id == IdToHuman::toTag(NonDimensional::PRANDTL))
      {
         return NonDimensional::PRANDTL;

      } else if(id == IdToHuman::toTag(NonDimensional::RAYLEIGH))
      {
         return NonDimensional::RAYLEIGH;

      } else if(id == IdToHuman::toTag(NonDimensional::ROBERTS))
      {
         return NonDimensional::ROBERTS;

      } else if(id == IdToHuman::toTag(NonDimensional::ROSSBY))
      {
         return NonDimensional::ROSSBY;

      } else if(id == IdToHuman::toTag(NonDimensional::RRATIO))
      {
         return NonDimensional::RRATIO;

      } else if(id == IdToHuman::toTag(NonDimensional::YXRATIO))
      {
         return NonDimensional::YXRATIO;

      } else if(id == IdToHuman::toTag(NonDimensional::ZXRATIO))
      {
         return NonDimensional::ZXRATIO;

      } else if(id == IdToHuman::toTag(NonDimensional::RO))
      {
         return NonDimensional::RO;

      } else if(id == IdToHuman::toTag(NonDimensional::TAYLOR))
      {
         return NonDimensional::TAYLOR;

      } else if(id == IdToHuman::toTag(NonDimensional::THETA))
      {
         return NonDimensional::THETA;

      } else
      {
         throw Exception("Unknown string to ID conversion requested (Nondimensional)");
      }
   }

   ModelOperator::Id HumanToId::toModOp(const std::string& id)
   {
      if(id == IdToHuman::toString(ModelOperator::QI))
      {
         return ModelOperator::QI;

      } else if(id == IdToHuman::toString(ModelOperator::TIME))
      {
         return ModelOperator::TIME;

      } else if(id == IdToHuman::toString(ModelOperator::IMPLICIT_LINEAR))
      {
         return ModelOperator::IMPLICIT_LINEAR;

      } else if(id == IdToHuman::toString(ModelOperator::EXPLICIT_LINEAR))
      {
         return ModelOperator::EXPLICIT_LINEAR;

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

      } else if(id == IdToHuman::toTag(PhysicalNames::DENSITY))
      {
         return PhysicalNames::DENSITY;

      } else if(id == IdToHuman::toTag(PhysicalNames::DZ_MEANTEMPERATURE))
      {
         return PhysicalNames::DZ_MEANTEMPERATURE;

      } else if(id == IdToHuman::toTag(PhysicalNames::ENTROPY))
      {
         return PhysicalNames::ENTROPY;

      } else if(id == IdToHuman::toTag(PhysicalNames::MAGNETIC))
      {
         return PhysicalNames::MAGNETIC;

      } else if(id == IdToHuman::toTag(PhysicalNames::MEANTEMPERATURE))
      {
         return PhysicalNames::MEANTEMPERATURE;

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

      } else if(id == IdToHuman::toTag(PhysicalNames::VORTICITY))
      {
         return PhysicalNames::VORTICITY;

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

      } else
      {
         throw Exception("Unknown string to ID conversion requested (PhysicalNames)");
      }
   }

   FieldComponents::Spectral::Id HumanToId::toComp(const std::string& id)
   {
      if(id == IdToHuman::toTag(FieldComponents::Spectral::ONE))
      {
         return FieldComponents::Spectral::ONE;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::TWO))
      {
         return FieldComponents::Spectral::TWO;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::THREE))
      {
         return FieldComponents::Spectral::THREE;
      } else if(id == IdToHuman::toTag(FieldComponents::Spectral::SCALAR))
      {
         return FieldComponents::Spectral::SCALAR;
      } else
      {
         throw Exception("Unknown string to ID conversion requested (FieldComponents::Spectral)");
      }
   }
}

}
