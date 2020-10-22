/**
 * @file ISphericalScalarLSpectrumWriter.cpp
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
#include "IoVariable/ISphericalScalarLSpectrumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/SpectrumTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalScalarLSpectrumWriter::ISphericalScalarLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalScalarEnergyBaseWriter(prefix + SpectrumTags::LBASENAME, SpectrumTags::EXTENSION, prefix + SpectrumTags::HEADER, type, SpectrumTags::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mEnergy(0)
   {
   }

   ISphericalScalarLSpectrumWriter::~ISphericalScalarLSpectrumWriter()
   {
   }

   void ISphericalScalarLSpectrumWriter::init()
   {
      this->mEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));

      ISphericalScalarEnergyBaseWriter::init();
   }

   void ISphericalScalarLSpectrumWriter::initializeEnergy()
   {
      this->mEnergy.setZero();
   }

   void ISphericalScalarLSpectrumWriter::storeEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mEnergy(l) += energy;
   }

   void ISphericalScalarLSpectrumWriter::write()
   {
      // Normalize by the volume
      this->mEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy.data(), this->mEnergy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << "#Time: "<< std::setprecision(14) << this->mTime << std::endl;

         // Total
         for(int i = 0; i < this->mEnergy.size(); i++)
         {
            this->mFile << i << "\t" << std::setprecision(14) << this->mEnergy(i) << std::endl;
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
