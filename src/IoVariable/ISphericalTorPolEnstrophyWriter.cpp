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
#include "IoVariable/EnstrophyTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolEnstrophyWriter::ISphericalTorPolEnstrophyWriter(const std::string& prefix, const std::string& type)
      : ISphericalTorPolEnstrophyBaseWriter(prefix + EnstrophyTags::BASENAME, EnstrophyTags::EXTENSION, prefix + EnstrophyTags::HEADER, type, EnstrophyTags::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mTorEnstrophy(2), mPolEnstrophy(2)
   {
      this->mTorEnstrophy.setConstant(-1);
      this->mPolEnstrophy.setConstant(-1);
   }

   ISphericalTorPolEnstrophyWriter::~ISphericalTorPolEnstrophyWriter()
   {
   }

   void ISphericalTorPolEnstrophyWriter::initializeEnstrophy()
   {
      this->mTorEnstrophy.setZero();
      this->mPolEnstrophy.setZero();
   }

   void ISphericalTorPolEnstrophyWriter::storePEnstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 0)
      {
         this->mPolEnstrophy(0) += enstrophy;
      } else
      {
         this->mPolEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeT1Enstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 0)
      {
         this->mTorEnstrophy(0) += enstrophy;
      } else
      {
         this->mTorEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::storeT2Enstrophy(const int l, const int m, const MHDFloat enstrophy)
   {
      if((l - m)%2 == 1)
      {
         this->mTorEnstrophy(0) += enstrophy;
      } else
      {
         this->mTorEnstrophy(1) += enstrophy;
      }
   }

   void ISphericalTorPolEnstrophyWriter::write()
   {
      // Normalize by the volume
      this->mTorEnstrophy /= this->mVolume;
      this->mPolEnstrophy /= this->mVolume;

      // Create file
      this->preWrite();

      // Get the "global" Kinetic enstrophy from MPI code
      #ifdef QUICC_MPI
         Array enstrophy(4);

         enstrophy(0) = this->mTorEnstrophy(0);
         enstrophy(1) = this->mTorEnstrophy(1);
         enstrophy(2) = this->mPolEnstrophy(0);
         enstrophy(3) = this->mPolEnstrophy(1);

         MPI_Allreduce(MPI_IN_PLACE, enstrophy.data(), enstrophy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnstrophy(0) = enstrophy(0);
         this->mTorEnstrophy(1) = enstrophy(1);
         this->mPolEnstrophy(0) = enstrophy(2);
         this->mPolEnstrophy(1) = enstrophy(3);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Total
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnstrophy.sum() + this->mPolEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mTorEnstrophy(0) + this->mPolEnstrophy(0) << "\t" << this->mTorEnstrophy(1) + this->mPolEnstrophy(1);
         }

         // Toroidal
         this->mFile << std::setprecision(14) << "\t" << this->mTorEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mTorEnstrophy(0) << "\t" << this->mTorEnstrophy(1);
         }

         // Poloidal
         this->mFile << std::setprecision(14) << "\t" << this->mPolEnstrophy.sum();
         if(this->mShowParity)
         {
            this->mFile << std::setprecision(14) << "\t" << this->mPolEnstrophy(0) << "\t" << this->mPolEnstrophy(1);
         }

         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic enstrophy is NaN
      if(std::isnan(this->mTorEnstrophy.sum()) || std::isnan(this->mPolEnstrophy.sum()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Toroidal/Poloidal enstrophy is NaN!");
      }
   }

}
}
