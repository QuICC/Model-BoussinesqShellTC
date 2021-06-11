/** 
 * @file ISphericalTorPolEnstrophyMSpectrumWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy m-spectrum calculation for a vector field in a sphere
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
#include "IoVariable/ISphericalTorPolEnstrophyMSpectrumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/EnstrophySpectrumTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolEnstrophyMSpectrumWriter::ISphericalTorPolEnstrophyMSpectrumWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyBaseWriter(prefix + EnstrophySpectrumTags::MBASENAME, EnstrophySpectrumTags::EXTENSION, prefix + EnstrophySpectrumTags::HEADER, type, EnstrophySpectrumTags::VERSION, Dimensions::Space::SPECTRAL, OVERWRITE), mTorEnstrophy(0), mPolEnstrophy(0)
   {
   }

   ISphericalTorPolEnstrophyMSpectrumWriter::~ISphericalTorPolEnstrophyMSpectrumWriter()
   {
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::init()
   {
      // Resize storage for spectra
      this->mTorEnstrophy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL));
      this->mPolEnstrophy = Array::Zero(this->res().sim()->dim(Dimensions::Simulation::SIM3D,Dimensions::Space::SPECTRAL));

      ISphericalTorPolEnstrophyBaseWriter::init();
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::initializeEnstrophy()
   {
      this->mTorEnstrophy.setZero();
      this->mPolEnstrophy.setZero();
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::storePEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mPolEnstrophy(l) += enstrophy; 
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::storeT1Enstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mTorEnstrophy(l) += enstrophy; 
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::storeT2Enstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      this->mTorEnstrophy(l) += enstrophy; 
   }

   void ISphericalTorPolEnstrophyMSpectrumWriter::write()
   {
      // Normalize by the volume
      this->mTorEnstrophy /= this->mVolume;
      this->mPolEnstrophy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic enstrophy from MPI code
      #ifdef QUICC_MPI
         Array enstrophy(2*this->mTorEnstrophy.size());

         enstrophy.segment(0,this->mTorEnstrophy.size()) = this->mTorEnstrophy;
         enstrophy.segment(this->mTorEnstrophy.size(),this->mPolEnstrophy.size()) = this->mPolEnstrophy;

         MPI_Allreduce(MPI_IN_PLACE, enstrophy.data(), enstrophy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnstrophy = enstrophy.segment(0,this->mTorEnstrophy.size());
         this->mPolEnstrophy = enstrophy.segment(this->mTorEnstrophy.size(),this->mPolEnstrophy.size());
      #endif //QUICC_MPI

      int ioPrec = 14;
      int ioIW = 5;
      int ioW = ioPrec+5;

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << "# time: " << std::setprecision(ioPrec);
         this->mFile << std::setw(ioW) << this->mTime << std::endl;

         this->mFile << "# enstrophy: " << std::setprecision(ioPrec);
         this->mFile << std::setw(ioW) << this->mTorEnstrophy.sum() + this->mPolEnstrophy.sum() << "\t";
         this->mFile << std::setw(ioW) << this->mTorEnstrophy.sum() << "\t";
         this->mFile << std::setw(ioW) << this->mPolEnstrophy.sum() << std::endl;

         this->mFile << std::left << std::setw(ioIW) << "# m" << "\t";
         this->mFile << std::setw(ioW) << "total" << "\t";
         this->mFile << std::setw(ioW) << "toroidal" << "\t";
         this->mFile << std::setw(ioW) << "poloidal" << std::endl;

         // Total
         for(int i = 0; i < this->mTorEnstrophy.size(); i++)
         {
            this->mFile << std::left << std::setw(ioIW) << i << "\t" << std::setprecision(ioPrec);
            this->mFile << std::setw(ioW) << this->mTorEnstrophy(i) + this->mPolEnstrophy(i) << "\t";
            this->mFile << std::setw(ioW) << this->mTorEnstrophy(i) << "\t";
            this->mFile << std::setw(ioW) << this->mPolEnstrophy(i) << std::endl;
         }
      }

      // Close file
      this->postWrite();

      // Abort if kinetic enstrophy is NaN
      if(std::isnan(this->mTorEnstrophy.sum()) || std::isnan(this->mPolEnstrophy.sum()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Toroidal/Poloidal L enstrophy spectrum is NaN!");
      }
   }

}
}
