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

         case PhysicalNames::DX_MEANTEMPERATURE:
            return "D_x mean temperature";

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

         case PhysicalNames::TILTED_TEMPERATURE:
            return "Tilted temperature";

         case PhysicalNames::TILTED_NO_STREAMFUNCTION:
            return "Tilted non orthogonal streamfunction";

         case PhysicalNames::TILTED_NO_VELOCITYZ:
            return "Tilted non orthogonal velocityZ";

         case PhysicalNames::TILTED_NO_VORTICITYZ:
            return "Tilted non orthogonal vorticityZ";

         case PhysicalNames::TILTED_STREAMFUNCTION:
            return "Tilted vertical streamfunction";

         case PhysicalNames::TILTED_VELOCITYZ:
            return "Tilted vertical velocityZ";

         case PhysicalNames::KINETIC_ENERGY:
            return "Kinetic energy";

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

         case PhysicalNames::DX_MEANTEMPERATURE:
            return "dx_meantemperature";

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

         case PhysicalNames::TILTED_TEMPERATURE:
            return "tilted_temperature";

         case PhysicalNames::TILTED_NO_STREAMFUNCTION:
            return "tilted_no_streamfunction";

         case PhysicalNames::TILTED_NO_VELOCITYZ:
            return "tilted_no_velocityz";

         case PhysicalNames::TILTED_NO_VORTICITYZ:
            return "tilted_no_vorticityz";

         case PhysicalNames::TILTED_STREAMFUNCTION:
            return "tilted_streamfunction";

         case PhysicalNames::TILTED_VELOCITYZ:
            return "tilted_velocityz";

         case PhysicalNames::KINETIC_ENERGY:
            return "kinetic_energy";

         default:
            throw Exception("Unknown ID to tag conversion requested (PhysicalNames)");
      }
   }

   std::string IdToHuman::toString(const FieldComponents::Physical::Id id)
   {
      switch(id)
      {
         case FieldComponents::Physical::X:
            return "X";
         case FieldComponents::Physical::Y:
            return "Y";
         case FieldComponents::Physical::Z:
            return "Z";
         case FieldComponents::Physical::R:
            return "R";
         case FieldComponents::Physical::THETA:
            return "Theta";
         case FieldComponents::Physical::PHI:
            return "Phi";
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
         case FieldComponents::Physical::X:
            return "x";
         case FieldComponents::Physical::Y:
            return "y";
         case FieldComponents::Physical::Z:
            return "z";
         case FieldComponents::Physical::R:
            return "r";
         case FieldComponents::Physical::THETA:
            return "theta";
         case FieldComponents::Physical::PHI:
            return "phi";
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
         case FieldComponents::Spectral::X:
            return "X";
         case FieldComponents::Spectral::Y:
            return "Y";
         case FieldComponents::Spectral::Z:
            return "Z";
         case FieldComponents::Spectral::R:
            return "R";
         case FieldComponents::Spectral::THETA:
            return "Theta";
         case FieldComponents::Spectral::PHI:
            return "Phi";
         case FieldComponents::Spectral::TOR:
            return "Toroidal";
         case FieldComponents::Spectral::POL:
            return "Poloidal";
         case FieldComponents::Spectral::Q:
            return "Q";
         case FieldComponents::Spectral::S:
            return "S";
         case FieldComponents::Spectral::T:
            return "T";
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
         case FieldComponents::Spectral::X:
            return "x";
         case FieldComponents::Spectral::Y:
            return "y";
         case FieldComponents::Spectral::Z:
            return "z";
         case FieldComponents::Spectral::R:
            return "r";
         case FieldComponents::Spectral::THETA:
            return "theta";
         case FieldComponents::Spectral::PHI:
            return "phi";
         case FieldComponents::Spectral::TOR:
            return "tor";
         case FieldComponents::Spectral::POL:
            return "pol";
         case FieldComponents::Spectral::Q:
            return "q";
         case FieldComponents::Spectral::S:
            return "s";
         case FieldComponents::Spectral::T:
            return "t";
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

         case NonDimensional::RO:
            return "ro";

         case NonDimensional::RRATIO:
            return "rratio";

         case NonDimensional::SCALE1D:
            return "scale1d";

         case NonDimensional::SCALE2D:
            return "scale2d";

         case NonDimensional::SCALE3D:
            return "scale3d";

         case NonDimensional::TAYLOR:
            return "taylor";

         case NonDimensional::THETA:
            return "theta";

         case NonDimensional::HEATING:
            return "heating";

         default:
            throw Exception("Unknown ID to tag conversion requested (NonDimensional)");
      }
   }

   std::string IdToHuman::toString(const ModelOperator::Id id)
   {
      switch(id)
      {
         case ModelOperator::TIME:
            return "time";

         case ModelOperator::IMPLICIT_LINEAR:
            return "implicit_linear";

         case ModelOperator::EXPLICIT_LINEAR:
            return "explicit_linear";

         case ModelOperator::EXPLICIT_NONLINEAR:
            return "explicit_nonlinear";

         case ModelOperator::EXPLICIT_NEXTSTEP:
            return "explicit_nextstep";

         case ModelOperator::STENCIL:
            return "stencil";

         default:
            throw Exception("Unknown ID to string conversion requested (ModelOperator)");
      }
   }

}
}
