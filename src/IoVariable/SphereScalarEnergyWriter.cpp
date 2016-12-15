/** 
 * @file SphereScalarEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field in a sphere
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
#include "IoVariable/SphereScalarEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   SphereScalarEnergyWriter::SphereScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-1.0)
   {
   }

   SphereScalarEnergyWriter::~SphereScalarEnergyWriter()
   {
   }

   void SphereScalarEnergyWriter::init()
   {
      // Normalize by sphere volume: 4/3*pi*r_o^3
      this->mVolume = (4.0/3.0)*Math::PI;

      IVariableAsciiEWriter::init();
   }

   void SphereScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);

      // Dealias variable data
      coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());
      
      // Recover dealiased BWD data
      Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

      // Get FWD storage
      Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

      // Compute projection transform for first dimension 
      coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

      // Compute |f|^2
      rOutVar.rData() = rOutVar.rData().array()*rOutVar.rData().conjugate().array();

      // Compute projection transform for first dimension 
      coord.transform1D().integrate_full(rInVar.rData(), rOutVar.data(), Transform::TransformCoordinatorType::Transform1DType::IntegratorType::ENERGYR2);

      // Compute integral over Chebyshev expansion and sum harmonics
      this->mEnergy = 0.0;

      #if defined GEOMHDISCC_SPATIALSCHEME_BLFM || defined GEOMHDISCC_SPATIALSCHEME_WLFM
         int start = 0;
         int  m0 = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0);
         if(m0 == 0)
         {
            this->mEnergy += (rInVar.slice(0).row(0).real()).sum();
            start = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT2D>(0);
         }

         this->mEnergy += 2.0*(rInVar.data().rightCols(rInVar.data().cols()-start).row(0).real()).sum();
      #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFM || defined GEOMHDISCC_SPATIALSCHEME_WLFM
      #if defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_WLFL
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int start = 0;
            if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,k) == 0)
            {
               this->mEnergy += (rInVar.slice(k)(0,0).real());
               start = 1;
            }
            this->mEnergy += 2.0*(rInVar.slice(k).rightCols(rInVar.slice(k).cols()-start).row(0).real()).sum();
         }
      #endif //defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_WLFL

      // Free BWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

      // Free FWD storage
      coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

      // Normalize by the spherical volume
      this->mEnergy /= this->mVolume;
   }

   void SphereScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &this->mEnergy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

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

         throw Exception("Scalar energy is NaN!");
      }
   }

}
}
