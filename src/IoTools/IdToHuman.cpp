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

         case PhysicalNames::ZONAL_VELOCITY:
            return "Zonal velocity";

         case PhysicalNames::NONZONAL_VELOCITY:
            return "Non zonal velocity";

         case PhysicalNames::VORTICITY:
            return "Vorticity";

         case PhysicalNames::VORTICITYX:
            return "VorticityX";

         case PhysicalNames::VORTICITYY:
            return "VorticityY";

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

         case PhysicalNames::FLUCT_TEMPERATURE:
            return "Fluctuating temperature";

         case PhysicalNames::FLUCT_MAGNETIC:
            return "Fluctuating magnetic";

         case PhysicalNames::FLUCT_MAGNETICX:
            return "Fluctuating magneticX";

         case PhysicalNames::FLUCT_MAGNETICY:
            return "Fluctuating magneticY";

         case PhysicalNames::FLUCT_MAGNETICZ:
            return "Fluctuating magneticZ";

         case PhysicalNames::FLUCT_VELOCITY:
            return "Fluctuating velocity";

         case PhysicalNames::FLUCT_VELOCITYX:
            return "Fluctuating velocityX";

         case PhysicalNames::FLUCT_VELOCITYY:
            return "Fluctuating velocityY";

         case PhysicalNames::FLUCT_VELOCITYZ:
            return "Fluctuating velocityZ";

         case PhysicalNames::MEAN_TEMPERATURE:
            return "Mean temperature";

         case PhysicalNames::MEAN_MAGNETIC:
            return "Mean magnetic";

         case PhysicalNames::MEAN_MAGNETICX:
            return "Mean magneticX";

         case PhysicalNames::MEAN_MAGNETICY:
            return "Mean magneticY";

         case PhysicalNames::MEAN_MAGNETICZ:
            return "Mean magneticZ";

         case PhysicalNames::MEAN_VELOCITY:
            return "Mean velocity";

         case PhysicalNames::MEAN_VELOCITYX:
            return "Mean velocityX";

         case PhysicalNames::MEAN_VELOCITYY:
            return "Mean velocityY";

         case PhysicalNames::MEAN_VELOCITYZ:
            return "Mean velocityZ";

         case PhysicalNames::KINETIC_ENERGY:
            return "Kinetic energy";

         case PhysicalNames::ZONAL_KINETIC_ENERGY:
            return "Zonal kinetic energy";

         case PhysicalNames::NONZONAL_KINETIC_ENERGY:
            return "Non zonal kinetic energy";

         case PhysicalNames::IMPOSED_MAGNETIC:
            return "Imposed magnetic field";

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

         case PhysicalNames::ZONAL_VELOCITY:
            return "zonal_velocity";

         case PhysicalNames::NONZONAL_VELOCITY:
            return "nonzonal_velocity";

         case PhysicalNames::VORTICITY:
            return "vorticity";

         case PhysicalNames::VORTICITYX:
            return "vorticityx";

         case PhysicalNames::VORTICITYY:
            return "vorticityy";

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

         case PhysicalNames::FLUCT_TEMPERATURE:
            return "fluct_temperature";

         case PhysicalNames::FLUCT_MAGNETIC:
            return "fluct_magnetic";

         case PhysicalNames::FLUCT_MAGNETICX:
            return "fluct_magneticx";

         case PhysicalNames::FLUCT_MAGNETICY:
            return "fluct_magneticy";

         case PhysicalNames::FLUCT_MAGNETICZ:
            return "fluct_magneticz";

         case PhysicalNames::FLUCT_VELOCITY:
            return "fluct_velocity";

         case PhysicalNames::FLUCT_VELOCITYX:
            return "fluct_velocityx";

         case PhysicalNames::FLUCT_VELOCITYY:
            return "fluct_velocityy";

         case PhysicalNames::FLUCT_VELOCITYZ:
            return "fluct_velocityz";

         case PhysicalNames::MEAN_TEMPERATURE:
            return "mean_temperature";

         case PhysicalNames::MEAN_MAGNETIC:
            return "mean_magnetic";

         case PhysicalNames::MEAN_MAGNETICX:
            return "mean_magneticx";

         case PhysicalNames::MEAN_MAGNETICY:
            return "mean_magneticy";

         case PhysicalNames::MEAN_MAGNETICZ:
            return "mean_magneticz";

         case PhysicalNames::MEAN_VELOCITY:
            return "mean_velocity";

         case PhysicalNames::MEAN_VELOCITYX:
            return "mean_velocityx";

         case PhysicalNames::MEAN_VELOCITYY:
            return "mean_velocityy";

         case PhysicalNames::MEAN_VELOCITYZ:
            return "mean_velocityz";

         case PhysicalNames::KINETIC_ENERGY:
            return "kinetic_energy";

         case PhysicalNames::ZONAL_KINETIC_ENERGY:
            return "zonal_kinetic_energy";

         case PhysicalNames::NONZONAL_KINETIC_ENERGY:
            return "nonzonal_kinetic_energy";

         case PhysicalNames::IMPOSED_MAGNETIC:
            return "imposed_magnetic";

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
         //
         // Nondimensional numbers
         case NonDimensional::EKMAN:
            return "ekman";

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

         case NonDimensional::TAYLOR:
            return "taylor";
        
         //
         // Geometrical numbers
         case NonDimensional::RO:
            return "ro";

         case NonDimensional::RRATIO:
            return "rratio";

         //
         // Flags
         case NonDimensional::HEATING:
            return "heating";

         //
         // Axis scalings
         case NonDimensional::SCALE1D:
            return "scale1d";

         case NonDimensional::SCALE2D:
            return "scale2d";

         case NonDimensional::SCALE3D:
            return "scale3d";

         // 
         // Greek alphabet
         case NonDimensional::ALPHA:
            return "alpha";

         case NonDimensional::BETA:
            return "beta";

         case NonDimensional::GAMMA:
            return "gamma";

         case NonDimensional::DELTA:
            return "delta";

         case NonDimensional::EPSILON:
            return "epsilon";

         case NonDimensional::ZETA:
            return "zeta";

         case NonDimensional::ETA:
            return "eta";

         case NonDimensional::THETA:
            return "theta";

         case NonDimensional::IOTA:
            return "iota";

         case NonDimensional::KAPPA:
            return "kappa";

         case NonDimensional::LAMBDA:
            return "lambda";

         case NonDimensional::MU:
            return "mu";

         case NonDimensional::NU:
            return "nu";

         case NonDimensional::XI:
            return "xi";

         case NonDimensional::OMICRON:
            return "omicron";

         case NonDimensional::PI:
            return "pi";

         case NonDimensional::RHO:
            return "rho";

         case NonDimensional::SIGMA:
            return "sigma";

         case NonDimensional::TAU:
            return "tau";

         case NonDimensional::UPSILON:
            return "upsilon";

         case NonDimensional::PHI:
            return "phi";

         case NonDimensional::CHI:
            return "chi";

         case NonDimensional::PSI:
            return "psi";

         case NonDimensional::OMEGA:
            return "omega";

         // 
         // Special flags
         case NonDimensional::ELEVATOR:
            return "elevator";

         case NonDimensional::FAST_MEAN:
            return "fast_mean";

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

         case ModelOperator::BOUNDARY:
            return "boundary";

         case ModelOperator::EXPLICIT_LINEAR:
            return "explicit_linear";

         case ModelOperator::EXPLICIT_NONLINEAR:
            return "explicit_nonlinear";

         case ModelOperator::EXPLICIT_NEXTSTEP:
            return "explicit_nextstep";

         case ModelOperator::STENCIL:
            return "stencil";

         case ModelOperator::INHOMOGENEOUS:
            return "inhomogeneous";

         default:
            throw Exception("Unknown ID to string conversion requested (ModelOperator)");
      }
   }

}
}
