/** 
 * @file ShellTorPolTorqueWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>

// External includes
//

// Class include
//
#include "IoVariable/ShellTorPolTorqueWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace IoVariable {

   ShellTorPolTorqueWriter::ShellTorPolTorqueWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   ShellTorPolTorqueWriter::~ShellTorPolTorqueWriter()
   {
   }

   void ShellTorPolTorqueWriter::init()
   {
      IVariableAsciiEWriter::init();
   }

   void ShellTorPolTorqueWriter::compute(Transform::TransformCoordinatorType& coord)
   {
	   // for the detection of the core carrying the Toroidal 1/0 function


   }

   void ShellTorPolTorqueWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(2);

         energy(0) = this->mTorEnergy;
         energy(1) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy(0);
         this->mPolEnergy = energy(1);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy) || std::isnan(this->mPolEnergy))
      {
         FrameworkMacro::abort(99);

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
