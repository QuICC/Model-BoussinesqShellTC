/**
 * @file ISphericalScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <iomanip>
#include <stdexcept>

// External includes
//

// Class include
//
#include "IoVariable/ISphericalScalarEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/EnergyTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalScalarEnergyWriter::ISphericalScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyBaseWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL, ISphericalScalarEnergyBaseWriter::EXTEND), mEnergy(2)
   {
      this->mEnergy.setConstant(-1);
   }

   ISphericalScalarEnergyWriter::~ISphericalScalarEnergyWriter()
   {
   }

   void ISphericalScalarEnergyWriter::initializeEnergy()
   {
      this->mEnergy.setZero();
   }

   void ISphericalScalarEnergyWriter::storeEnergy(const int l, const int m, const MHDFloat energy)
   {
      if((l - m)%2 == 0)
      {
         this->mEnergy(0) += energy;
      } else
      {
         this->mEnergy(1) += energy;
      }
   }

   void ISphericalScalarEnergyWriter::write()
   {
      // Normalize by the volume
      this->mEnergy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mEnergy(0) << "\t" << this->mEnergy(1);
         }
         
         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy.sum()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Scalar energy is NaN!");
      }
   }

}
}
