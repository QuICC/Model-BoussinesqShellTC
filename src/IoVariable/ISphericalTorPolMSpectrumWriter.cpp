/** 
 * @file ISphericalTorPolMSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
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
#include "IoVariable/ISphericalTorPolMSpectrumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/SpectrumTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolMSpectrumWriter::ISphericalTorPolMSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyBaseWriter(prefix + SpectrumTags::MBASENAME, SpectrumTags::EXTENSION, prefix + SpectrumTags::HEADER, type, SpectrumTags::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorEnergy(0), mPolEnergy(0)
   {
   }

   ISphericalTorPolMSpectrumWriter::~ISphericalTorPolMSpectrumWriter()
   {
   }

   void ISphericalTorPolMSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL));
      this->mPolEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolEnergyBaseWriter::init();
   }

   void ISphericalTorPolMSpectrumWriter::initializeEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ISphericalTorPolMSpectrumWriter::storeQEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(m) += energy; 
   }

   void ISphericalTorPolMSpectrumWriter::storeSEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(m) += energy; 
   }

   void ISphericalTorPolMSpectrumWriter::storeTEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mTorEnergy(m) += energy; 
   }

   void ISphericalTorPolMSpectrumWriter::write()
   {
      // Normalize by the volume
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(2*this->mTorEnergy.size());

         energy.segment(0,this->mTorEnergy.size()) = this->mTorEnergy;
         energy.segment(this->mTorEnergy.size(),this->mPolEnergy.size()) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy.segment(0,this->mTorEnergy.size());
         this->mPolEnergy = energy.segment(this->mTorEnergy.size(),this->mPolEnergy.size());
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
         this->mFile << std::setw(ioW) << this->mTorEnergy.sum() + this->mPolEnergy.sum() << "\t";
         this->mFile << std::setw(ioW) << this->mTorEnergy.sum() << "\t";
         this->mFile << std::setw(ioW) << this->mPolEnergy.sum() << std::endl;

         this->mFile << std::left << std::setw(ioIW) << "# m" << "\t";
         this->mFile << std::setw(ioW) << "total" << "\t";
         this->mFile << std::setw(ioW) << "toroidal" << "\t";
         this->mFile << std::setw(ioW) << "poloidal" << std::endl;

         // Total
         for(int i = 0; i < this->mTorEnergy.size(); i++)
         {
            this->mFile << std::left << std::setw(ioIW) << i << "\t" << std::setprecision(ioPrec);
            this->mFile << std::setw(ioW) << this->mTorEnergy(i) + this->mPolEnergy(i) << "\t";
            this->mFile << std::setw(ioW) << this->mTorEnergy(i) << "\t";
            this->mFile << std::setw(ioW) << this->mPolEnergy(i) << std::endl;
         }
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy.sum()) || std::isnan(this->mPolEnergy.sum()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Toroidal/Poloidal L energy spectrum is NaN!");
      }
   }

}
}
