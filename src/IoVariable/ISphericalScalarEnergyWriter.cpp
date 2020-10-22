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
      : IVariableAsciiWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL, IVariableAsciiWriter::EXTEND), mHasMOrdering(false), mEnergy(-1.0)
   {
   }

   ISphericalScalarEnergyWriter::~ISphericalScalarEnergyWriter()
   {
   }

   void ISphericalScalarEnergyWriter::init()
   {
      IVariableAsciiWriter::init();
   }

   void ISphericalScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Dealias variable data
      coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute projection transform for first dimension
      Array spectrum;
      coord.transform1D().integrate_energy(spectrum, rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::ENERGY_PROJ, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::ENERGY_R2);

      // Compute integral over Worland expansion and sum harmonics
      this->mEnergy = 0.0;

      if(this->mHasMOrdering)
      {
         int start = 0;
         int  m0 = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         if(m0 == 0)
         {
            int n = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0);
            this->mEnergy += spectrum.segment(0, n).sum();
            start += n;
         }

         this->mEnergy += 2.0*spectrum.segment(start, spectrum.size() - start).sum();
      } else
      {
         int start = 0;
         int n = 0;
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            n = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k);
            if(this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mEnergy += spectrum(start);
               start += 1;
               n -= 1;
            }
            this->mEnergy += 2.0*spectrum.segment(start,n).sum();
            start += n;
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Normalize by the spherical volume
      this->mEnergy /= this->mVolume;
   }

   void ISphericalScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mEnergy))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Scalar energy is NaN!");
      }
   }

}
}
