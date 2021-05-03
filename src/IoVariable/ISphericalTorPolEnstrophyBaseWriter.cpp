/**
 * @file ISphericalTorPolEnstrophyBaseWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics enstrophy calculation for a vector field in a spherical geometry
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
#include "IoVariable/ISphericalTorPolEnstrophyBaseWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/EnergyTags.hpp"

namespace QuICC {

namespace IoVariable {

   ISphericalTorPolEnstrophyBaseWriter::ISphericalTorPolEnstrophyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const ISphericalTorPolEnstrophyBaseWriter::WriteMode mode)
      : IVariableAsciiWriter(name, ext, header, type, version, id, mode), mHasMOrdering(false), mVolume(std::numeric_limits<MHDFloat>::quiet_NaN()), mShowParity(false)
   {
   }

   ISphericalTorPolEnstrophyBaseWriter::~ISphericalTorPolEnstrophyBaseWriter()
   {
   }

   void ISphericalTorPolEnstrophyBaseWriter::showParity()
   {
      this->mShowParity = true;
   }

   void ISphericalTorPolEnstrophyBaseWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator vIt;
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      Array spectrum;

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute first part of toroidal enstrophy integral for first dimension
      coord.transform1D().integrate_energy(spectrum, rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::ENERGY_PROJ, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::ENERGY_INTG);

      this->initializeEnergy();

      MHDFloat lfactor = 0.0;
      MHDFloat factor = 1.0;
      int idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = std::pow(l_*(l_+1.0),2);

               this->storeT1Enstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::pow(l_*(l_+1.0),2);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               // m = 0, no factor of two
               if(m_ == 0)
               {
                  factor = 1.0;
               } else
               {
                  factor = 2.0;
               }

               this->storeT1Enstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor2 = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute second part of toroidal enstrophy integral for first dimension
      coord.transform1D().integrate_energy(spectrum, rInVarTor2.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::ENERGY_DIFFR, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::ENERGY_INTG);

      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            // m = 0, no factor of two
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l_*(l_+1.0);

               this->storeT2Enstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l_*(l_+1.0);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               // m = 0, no factor of two
               if(m_ == 0)
               {
                  factor = 1.0;
               } else
               {
                  factor = 2.0;
               }

               this->storeT2Enstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
	}

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor2);

      // Dealias poloidal variable data 
      coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::POL));

      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarPolQ = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute poloidal enstrophy integral for first dimension 
      coord.transform1D().integrate_energy(spectrum, rInVarPolQ.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::ENERGY_SLAPL, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::ENERGY_R2);

      // Compute poloidal component of the enstrophy
      idx = 0;
      if(this->mHasMOrdering)
      {
         // Loop over harmonic order m
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            // m = 0, no factor of two
            int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            if(m_ == 0)
            {
               factor = 1.0;
            } else
            {
               factor = 2.0;
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = l_*(l_+1.0);

               this->storePEnstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = l_*(l_+1.0);
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               // m = 0, no factor of two
               if(m_ == 0)
               {
                  factor = 1.0;
               } else
               {
                  factor = 2.0;
               }

               this->storePEnstrophy(l_, m_, factor*lfactor*spectrum(idx));
               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarPolQ);
    }
}
}
