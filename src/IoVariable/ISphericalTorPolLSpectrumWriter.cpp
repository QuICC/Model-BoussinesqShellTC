/** 
 * @file ISphericalTorPolLSpectrumWriter.cpp
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
#include "IoVariable/ISphericalTorPolLSpectrumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/SpectrumTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolLSpectrumWriter::ISphericalTorPolLSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnergyBaseWriter(prefix + SpectrumTags::LBASENAME, SpectrumTags::EXTENSION, prefix + SpectrumTags::HEADER, type, SpectrumTags::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorEnergy(0), mPolEnergy(0)
   {
   }

   ISphericalTorPolLSpectrumWriter::~ISphericalTorPolLSpectrumWriter()
   {
   }

   void ISphericalTorPolLSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
      this->mPolEnergy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM2D,Dimensions::Space::SPECTRAL));
   }

   void ISphericalTorPolLSpectrumWriter::initializeEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ISphericalTorPolLSpectrumWriter::storeQEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(l) += energy; 
   }

   void ISphericalTorPolLSpectrumWriter::storeSEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mPolEnergy(l) += energy; 
   }

   void ISphericalTorPolLSpectrumWriter::storeTEnergy(const int l, const int m, const MHDFloat energy)
   {
      this->mTorEnergy(l) += energy; 
   }

   void ISphericalTorPolLSpectrumWriter::write()
   {
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

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << "#Time: "<< std::setprecision(14) << this->mTime << std::endl;

         // Total
         for(int i = 0; i < this->mTorEnergy.size(); i++)
         {
            this->mFile << i << "\t" << std::setprecision(14) << this->mTorEnergy(i) + this->mPolEnergy(i) << "\t" << this->mTorEnergy(i) << "\t" << this->mPolEnergy(i) << std::endl;
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
