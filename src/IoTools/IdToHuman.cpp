/** \file IdToHuman.cpp
 *  \brief Source of enum id to human strings converters
 */

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
         case PhysicalNames::PRESSURE:
            return "Pressure";
         case PhysicalNames::TEMPERATURE:
            return "Temperature";
         case PhysicalNames::STREAMFUNCTION:
            return "Streamfunction";
         case PhysicalNames::VELOCITYZ:
            return "VelocityZ";
         case PhysicalNames::VORTICITYZ:
            return "VorticityZ";
         case PhysicalNames::PHI:
            return "Phi";
         case PhysicalNames::VELOCITY:
            return "Velocity";
         case PhysicalNames::MAGNETIC:
            return "Magnetic";
         case PhysicalNames::VORTICITY:
            return "Vorticity";
         default:
            return "UNKNOWN";
      }
   }

   std::string IdToHuman::toTag(const PhysicalNames::Id id)
   {
      switch(id)
      {
         case PhysicalNames::CODENSITY:
            return "codensity";
         case PhysicalNames::PRESSURE:
            return "pressure";
         case PhysicalNames::TEMPERATURE:
            return "temperature";
         case PhysicalNames::STREAMFUNCTION:
            return "streamfunction";
         case PhysicalNames::VELOCITYZ:
            return "velocityz";
         case PhysicalNames::VORTICITYZ:
            return "vorticityz";
         case PhysicalNames::PHI:
            return "phi";
         case PhysicalNames::VELOCITY:
            return "velocity";
         case PhysicalNames::MAGNETIC:
            return "magnetic";
         case PhysicalNames::VORTICITY:
            return "vorticity";
         default:
            return "unknown";
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
            return "UNKNOWN";
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
            return "UNKNOWN";
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
            return "UNKNOWN";
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
            return "UNKNOWN";
      }
   }

}
}
