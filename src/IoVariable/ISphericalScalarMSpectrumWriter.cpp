/**
 * @file ISphericalScalarMSpectrumWriter.cpp
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
#include "IoVariable/ISphericalScalarMSpectrumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/SpectrumTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalScalarMSpectrumWriter::ISphericalScalarMSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyBaseWriter(prefix + SpectrumTags::MBASENAME, SpectrumTags::EXTENSION, prefix + SpectrumTags::HEADER, type, SpectrumTags::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mEnergy(0)
   {
   }

   ISphericalScalarMSpectrumWriter::~ISphericalScalarMSpectrumWriter()
   {
   }

   void ISphericalScalarMSpectrumWriter::init()
   {
      this->mEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL));

      ISphericalScalarEnergyBaseWriter::init();
   }

   void ISphericalScalarMSpectrumWriter::initializeEnergy()
   {
      this->mEnergy.setZero();
   }

   void ISphericalScalarMSpectrumWriter::storeEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mEnergy(m) += energy;
   }

   void ISphericalScalarMSpectrumWriter::write()
   {
      // Normalize by the volume
      this->mEnergy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      int ioPrec = 14;
      int ioIW = 5;
      int ioW = ioPrec+5;

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << "# time: " << std::setprecision(ioPrec);
         this->mFile << std::setw(ioW) << this->mTime << std::endl;

         this->mFile << "# energy: " << std::setprecision(ioPrec);
         this->mFile << std::setw(ioW) << this->mEnergy.sum() << std::endl;

         this->mFile << std::left << std::setw(ioIW) << "# m" << "\t";
         this->mFile << std::setw(ioW) << "total" << std::endl;

         // Total
         for(int i = 0; i < this->mEnergy.size(); i++)
         {
            this->mFile << std::left << std::setw(ioIW) << i << "\t" << std::setprecision(ioPrec);
            this->mFile << std::setw(ioW) << this->mEnergy(i) << std::endl;
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
