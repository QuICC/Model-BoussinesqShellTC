/**
 * @file SphereAngularMomentumWriter.cpp
 * @brief Source of the implementation of the ASCII sphere angular momentum
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
#include "IoVariable/SphereAngularMomentumWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/AngularMomentumTags.hpp"

namespace QuICC {

namespace IoVariable {

   SphereAngularMomentumWriter::SphereAngularMomentumWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiWriter(prefix + AngularMomentumTags::BASENAME, AngularMomentumTags::EXTENSION, prefix + AngularMomentumTags::HEADER, type, AngularMomentumTags::VERSION, Dimensions::Space::SPECTRAL, EXTEND), mMomentum(3)
   {
   }

   SphereAngularMomentumWriter::~SphereAngularMomentumWriter()
   {
   }

   void SphereAngularMomentumWriter::init()
   {
      #if defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
         this->mHasMOrdering = true;
      #endif //defined QUICC_SPATIALSCHEME_BLFM || defined QUICC_SPATIALSCHEME_WLFM
      #if defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL
         this->mHasMOrdering = false;
      #endif // defined QUICC_SPATIALSCHEME_BLFL || defined QUICC_SPATIALSCHEME_WLFL

      IVariableAsciiWriter::init();
   }

   void SphereAngularMomentumWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      // get iterator to field
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);

      ArrayZ mom;
      this->mMomentum.setZero();

      // Dealias toroidal variable data
      coord.communicator().dealiasSpectral(vRange.first->second->dom(0).total().comp(FieldComponents::Spectral::TOR));
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Compute energy integral for first dimension
      coord.transform1D().integrate_volume(mom, rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::VOLUME_PROJ, Transform::TransformCoordinatorType::Transform1DType::IntegratorType::VOLUME_R3);

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
               factor = std::sqrt(2.0);
            }

            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
               lfactor = std::sqrt(16.0*Math::PI/(2.0*l_+1.0));

               if(l_ == 1 && m_ == 0)
               {
                  this->mMomentum(2) = factor*lfactor*mom(idx).real();
               } else if(l_ == 1 && m_ == 1)
               {
                  this->mMomentum(0) = -factor*lfactor*mom(idx).real();
                  this->mMomentum(1) = factor*lfactor*mom(idx).imag();
               }

               idx += 1;
            }
         }
      } else
      {
         // Loop over harmonic degree l
         for(int k = 0; k < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int l_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
            lfactor = std::sqrt(16.0*Math::PI/(2.0*l_+1.0));
            for(int j = 0; j < this->res().cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
            {
               int m_ = this->res().cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j,k);
               if(l_ == 1 && m_ == 0)
               {
                  factor = 1.0;
                  this->mMomentum(2) = factor*lfactor*mom(idx).real();
               } else if(l_ == 1 && m_ == 1)
               {
                  factor = std::sqrt(2.0);
                  this->mMomentum(0) = -factor*lfactor*mom(idx).real();
                  this->mMomentum(1) = factor*lfactor*mom(idx).imag();
               }

               idx += 1;
            }
         }
      }

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);
   }

   void SphereAngularMomentumWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" value from MPI code
      #ifdef QUICC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mMomentum.data(), this->mMomentum.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mMomentum.norm() << "\t" << this->mMomentum.transpose();
         
         // End line
         this->mFile << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if is NaN
      if(std::isnan(this->mMomentum.norm()))
      {
         FrameworkMacro::abort(99);

         throw std::logic_error("Sphere angular momentum is NaN!");
      }
   }

}
}
