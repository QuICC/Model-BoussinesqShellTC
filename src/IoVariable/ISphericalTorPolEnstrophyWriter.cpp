/** 
 * @file ISphericalTorPolEnstrophyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy calculation for vector field in a spherical shell
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
#include "IoVariable/ISphericalTorPolEnstrophyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/EnergyTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolEnstrophyWriter::ISphericalTorPolEnstrophyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyBaseWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mTorEnergy(2), mPolEnergy(2)
   {
      this->mTorEnergy.setConstant(-1);
      this->mPolEnergy.setConstant(-1);
   }

   ISphericalTorPolEnstrophyWriter::~ISphericalTorPolEnstrophyWriter()
   {
   }

   void ISphericalTorPolEnstrophyWriter::initializeEnergy()
   {
      this->mTorEnergy.setZero();
      this->mPolEnergy.setZero();
   }

   void ISphericalTorPolEnstrophyWriter::storePEnstrophy(const int l, const int m, const MHDFloat energy)
   {
      if((l - m)%2 == 0)
      {
         this->mPolEnergy(0) += energy;
      } else
      {
         this->mPolEnergy(1) += energy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeT1Enstrophy(const int l, const int m, const MHDFloat energy)
   {
      if((l - m)%2 == 0)
      {
         this->mTorEnergy(0) += energy;
      } else
      {
         this->mTorEnergy(1) += energy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeT2Enstrophy(const int l, const int m, const MHDFloat energy)
   {
      if((l - m)%2 == 1)
      {
         this->mTorEnergy(0) += energy; 
      } else
      {
         this->mTorEnergy(1) += energy; 
      }
   }

   void ISphericalTorPolEnstrophyWriter::write()
   {
      // Normalize by the volume
      this->mTorEnergy /= 2.0*this->mVolume;
      this->mPolEnergy /= 2.0*this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(4);

         energy(0) = this->mTorEnergy(0);
         energy(1) = this->mTorEnergy(1);
         energy(2) = this->mPolEnergy(0);
         energy(3) = this->mPolEnergy(1);

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy(0) = energy(0);
         this->mTorEnergy(1) = energy(1);
         this->mPolEnergy(0) = energy(2);
         this->mPolEnergy(1) = energy(3);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Total 
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy.sum() + this->mPolEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mTorEnergy(0) + this->mPolEnergy(0) << "\t" << this->mTorEnergy(1) + this->mPolEnergy(1);
         }

         // Toroidal 
         this->mFile << std::setprecision(14) << "\t" << this->mTorEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mTorEnergy(0) << "\t" << this->mTorEnergy(1);
         }

         // Poloidal 
         this->mFile << std::setprecision(14) << "\t" << this->mPolEnergy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mPolEnergy(0) << "\t" << this->mPolEnergy(1);
         }
         
         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy.sum()) || std::isnan(this->mPolEnergy.sum()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
