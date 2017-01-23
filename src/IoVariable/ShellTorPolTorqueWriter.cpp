/** 
 * @file ShellTorPolTorqueWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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
#include "IoVariable/ShellTorPolTorqueWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoTools/IdToHuman.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace QuICC {

namespace IoVariable {

   ShellTorPolTorqueWriter::ShellTorPolTorqueWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mTorEnergy(-1.0), mPolEnergy(-1.0)
   {
   }

   ShellTorPolTorqueWriter::~ShellTorPolTorqueWriter()
   {
   }

   void ShellTorPolTorqueWriter::init()
   {
      IVariableAsciiEWriter::init();
   }

   void ShellTorPolTorqueWriter::compute(Transform::TransformCoordinatorType& coord)
   {

	   // for the detection of the core carrying the Toroidal 1/0 function


	   /*
	    * compute the "bad way"
	    */

	  // get iterator to field
	  vector_iterator vIt;
	  vector_iterator_range vRange = this->vectorRange();
	  assert(std::distance(vRange.first, vRange.second) == 1);
	  assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
	  assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

	  // Dealias toroidal variable data
	  coord.communicator().dealiasSpectral(vRange.first->second->rDom(0).rTotal().rComp(FieldComponents::Spectral::TOR));

	  // Recover dealiased BWD data
	  Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

	  // Get FWD storage
	  Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVarTor = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

	  // Compute projection transform for first dimension
	  coord.transform1D().project(rOutVarTor.rData(), rInVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

	  // Compute |f|^2
	  rOutVarTor.rData() = rOutVarTor.rData().array()*rOutVarTor.rData().conjugate().array();

	  // Compute projection transform for first dimension
	  coord.transform1D().integrate_full(rInVarTor.rData(), rOutVarTor.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::INTG);

	  // Compute integral over Chebyshev expansion and sum harmonics

	  MHDFloat lfactor = 0.0;
	  #ifdef QUICC_SPATIALSCHEME_SLFM
		 double factor = 1.0;
		 // Loop over harmonic order m
		 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		 {
			// determine current m
			int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
			// m = 0, no factor of two
			if( m == 0)
			{
			   factor = 1.0;
			} else
			{
			   factor = 2.0;
			}

			for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++)
			{
			   int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
			   lfactor = l*(l+1.0);

			   this->mTorEnergy(l,m) += factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
			   this->mTorRadial += factor*lfactor*(rOutVarTor.slice(k).col(j).real());
			}
		 }
	  #endif //defined QUICC_SPATIALSCHEME_SLFM
	  #ifdef QUICC_SPATIALSCHEME_SLFL
		 double factor=1.;
		 // Loop over harmonic degree l
		 for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
		 {
			int l = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(k);
			lfactor = l*(l+1.0);

			for(int j = 0; j < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(k); j++){
				int m = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(j, k);
				if(m==0){
					factor = 1.0;
				} else {
					factor = 2.0;
				}
				this->mTorEnergy(l,m) += factor*lfactor*(this->mSphIntgOp*rInVarTor.slice(k).col(j).real()).sum();
				this->mTorRadial += factor*lfactor*(rOutVarTor.slice(k).col(j).real());
			}

		 }
	  #endif //QUICC_SPATIALSCHEME_SLFL

	  // Free BWD storage
	  coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVarTor);

	  // Free FWD storage
	  coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVarTor);

	  // Normalize by sphere volume: 4/3*pi*(r_o^3 - r_i^3)
	  this->mTorEnergy /= 2*this->mVolume;

   }

   void ShellTorPolTorqueWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef QUICC_MPI
         Array energy(2);

         energy(0) = this->mTorEnergy;
         energy(1) = this->mPolEnergy;

         MPI_Allreduce(MPI_IN_PLACE, energy.data(), 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

         this->mTorEnergy = energy(0);
         this->mPolEnergy = energy(1);
      #endif //QUICC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mTorEnergy + this->mPolEnergy << "\t" << this->mTorEnergy << "\t" << this->mPolEnergy << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(this->mTorEnergy) || std::isnan(this->mPolEnergy))
      {
         FrameworkMacro::abort(99);

         throw Exception("Toroidal/Poloidal energy is NaN!");
      }
   }

}
}
