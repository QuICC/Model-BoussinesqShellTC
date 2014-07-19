/** 
 * @file IdToHuman.cpp
 * @brief Source of enum id to human strings converters
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
#include "IoTools/IdToHuman.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoTools {

   std::string IdToHuman::toString(const PhysicalNames::Id id)
   {
      switch(id)
      {
         case PhysicalNames::CODENSITY:
            return "Codensity";

         case PhysicalNames::DENSITY:
            return "Density";

         case PhysicalNames::DZ_MEANTEMPERATURE:
            return "D_z mean temperature";

         case PhysicalNames::ENTROPY:
            return "Entropy";

         case PhysicalNames::MAGNETIC:
            return "Magnetic";

         case PhysicalNames::MEANTEMPERATURE:
            return "Mean temperature";

         case PhysicalNames::PRESSURE:
            return "Pressure";

         case PhysicalNames::STREAMFUNCTION:
            return "Streamfunction";

         case PhysicalNames::TEMPERATURE:
            return "Temperature";

         case PhysicalNames::VELOCITY:
            return "Velocity";

         case PhysicalNames::VELOCITYX:
            return "VelocityX";

         case PhysicalNames::VELOCITYY:
            return "VelocityY";

         case PhysicalNames::VELOCITYZ:
            return "VelocityZ";

         case PhysicalNames::VORTICITY:
            return "Vorticity";

         case PhysicalNames::VORTICITYZ:
            return "VorticityZ";

         case PhysicalNames::PHI:
            return "Phi";

         case PhysicalNames::NO_STREAMFUNCTION:
            return "Non orthogonal streamfunction";

         case PhysicalNames::NO_VELOCITYZ:
            return "Non orthogonal velocityZ";

         case PhysicalNames::NO_VORTICITYZ:
            return "Non orthogonal vorticityZ";

         default:
            throw Exception("Unknown ID to string conversion requested (PhysicalNames)");
      }
   }

   std::string IdToHuman::toTag(const PhysicalNames::Id id)
   {
      switch(id)
      {
         case PhysicalNames::CODENSITY:
            return "codensity";

         case PhysicalNames::DENSITY:
            return "density";

         case PhysicalNames::DZ_MEANTEMPERATURE:
            return "dz_meantemperature";

         case PhysicalNames::ENTROPY:
            return "entropy";

         case PhysicalNames::MAGNETIC:
            return "magnetic";

         case PhysicalNames::MEANTEMPERATURE:
            return "meantemperature";

         case PhysicalNames::PRESSURE:
            return "pressure";

         case PhysicalNames::STREAMFUNCTION:
            return "streamfunction";

         case PhysicalNames::TEMPERATURE:
            return "temperature";

         case PhysicalNames::VELOCITY:
            return "velocity";

         case PhysicalNames::VELOCITYX:
            return "velocityx";

         case PhysicalNames::VELOCITYY:
            return "velocityy";

         case PhysicalNames::VELOCITYZ:
            return "velocityz";

         case PhysicalNames::VORTICITY:
            return "vorticity";

         case PhysicalNames::VORTICITYZ:
            return "vorticityz";

         case PhysicalNames::PHI:
            return "phi";

         case PhysicalNames::NO_STREAMFUNCTION:
            return "no_streamfunction";

         case PhysicalNames::NO_VELOCITYZ:
            return "no_velocityz";

         case PhysicalNames::NO_VORTICITYZ:
            return "no_vorticityz";

         default:
            throw Exception("Unknown ID to tag conversion requested (PhysicalNames)");
      }
   }

   std::string IdToHuman::toString(const FieldComponents::Physical::Id id)
   {
      switch(id)
      {
         case FieldComponents::Physical::ONE:
            return "X";
         case FieldComponents::Physical::TWO:
            return "Y";
         case FieldComponents::Physical::THREE:
            return "Z";
         case FieldComponents::Physical::SCALAR:
            return "";
         default:
            throw Exception("Unknown ID to string conversion requested (FieldComponents::Physical)");
      }
   }

   std::string IdToHuman::toTag(const FieldComponents::Physical::Id id)
   {
      switch(id)
      {
         case FieldComponents::Physical::ONE:
            return "x";
         case FieldComponents::Physical::TWO:
            return "y";
         case FieldComponents::Physical::THREE:
            return "z";
         case FieldComponents::Physical::SCALAR:
            return "";
         default:
            throw Exception("Unknown ID to tag conversion requested (FieldComponents::Physical)");
      }
   }

   std::string IdToHuman::toString(const FieldComponents::Spectral::Id id)
   {
      switch(id)
      {
         case FieldComponents::Spectral::ONE:
            return "Toroidal";
         case FieldComponents::Spectral::TWO:
            return "Poloidal";
         case FieldComponents::Spectral::THREE:
            return "NOTYET";
         case FieldComponents::Spectral::SCALAR:
            return "";
         default:
            throw Exception("Unknown ID to string conversion requested (FieldComponents::Spectral)");
      }
   }

   std::string IdToHuman::toTag(const FieldComponents::Spectral::Id id)
   {
      switch(id)
      {
         case FieldComponents::Spectral::ONE:
            return "tor";
         case FieldComponents::Spectral::TWO:
            return "pol";
         case FieldComponents::Spectral::THREE:
            return "NOTYET";
         case FieldComponents::Spectral::SCALAR:
            return "";
         default:
            throw Exception("Unknown ID to tag conversion requested (FieldComponents::Spectral)");
      }
   }

   std::string IdToHuman::toTag(const NonDimensional::Id id)
   {
      switch(id)
      {
         case NonDimensional::CHI:
            return "chi";

         case NonDimensional::EKMAN:
            return "ekman";

         case NonDimensional::GAMMA:
            return "gamma";

         case NonDimensional::GAPWIDTH:
            return "gapwidth";

         case NonDimensional::MAGEKMAN:
            return "magnetic_ekman";

         case NonDimensional::MAGPRANDTL:
            return "magnetic_prandtl";

         case NonDimensional::PRANDTL:
            return "prandtl";

         case NonDimensional::RAYLEIGH:
            return "rayleigh";

         case NonDimensional::ROBERTS:
            return "roberts";

         case NonDimensional::ROSSBY:
            return "rossby";

         case NonDimensional::RRATIO:
            return "rratio";

         case NonDimensional::TAYLOR:
            return "taylor";

         case NonDimensional::THETA:
            return "theta";

         default:
            throw Exception("Unknown ID to tag conversion requested (NonDimensional)");
      }
   }

   std::string IdToHuman::toString(const ModelOperator::Id id)
   {
      switch(id)
      {
         case ModelOperator::QI:
            return "qi";

         case ModelOperator::TIME:
            return "time";

         case ModelOperator::IMPLICIT_LINEAR:
            return "implicit_linear";

         case ModelOperator::EXPLICIT_LINEAR:
            return "explicit_linear";

         case ModelOperator::STENCIL:
            return "stencil";

         default:
            throw Exception("Unknown ID to string conversion requested (ModelOperator)");
      }
   }

}
}
