/** 
 * @file Cartesian1DScalarKurtWriter.cpp
 * @brief Source of the implementation of the ASCII Cartesian 1D (double periodic) statistics RMS calculation for scalar field
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
#include "IoStats/Cartesian1DScalarKurtWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

   namespace IoStats {

      Cartesian1DScalarKurtWriter::Cartesian1DScalarKurtWriter(const std::string& prefix, const std::string& type)
         : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-Array::Ones(2))
      {
      }

      Cartesian1DScalarKurtWriter::~Cartesian1DScalarKurtWriter()
      {
      }

      void Cartesian1DScalarKurtWriter::init()
      {

         // Here is where we get the mean from the 0th mode of the variable in spectral space 

         int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

         // Normalize by Cartesian Area A = Nx*Ny (need to worry about 2 pi in fft?)
         this->mArea = cols*cols;

         IVariableAsciiEWriter::init();
      }

      void Cartesian1DScalarKurtWriter::precompute(Transform::TransformCoordinatorType& coord)
      {
         //compute the mean here
         this->mAvg.setConstant(0.0);
         // Dealias variable data
         coord.communicator().dealiasSpectral(sRange.first->second->rDom(0).rTotal());

         // Recover dealiased BWD data
         Transform::TransformCoordinatorType::CommunicatorType::Bwd1DType &rInVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().recoverBwd();

         // Get FWD storage
         Transform::TransformCoordinatorType::CommunicatorType::Fwd1DType &rOutVar = coord.communicator().storage<Dimensions::Transform::TRA1D>().provideFwd();

         // Compute projection transform for first dimension 
         coord.transform1D().project(rOutVar.rData(), rInVar.data(), Transform::TransformCoordinatorType::Transform1DType::ProjectorType::PROJ);

         // set mAvg to be the 0th mode of the field
         this->mAvg = rOutVar.slice(0).col(0).real();

         // Free BWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeBwd(rInVar);

         // Free FWD storage
         coord.communicator().storage<Dimensions::Transform::TRA1D>().freeFwd(rOutVar);

      }

      void Cartesian1DScalarKurtWriter::compute(Transform::TransformCoordinatorType& coord)
      {
         // calculate the transforms and calculate the RMS

         // Initialize the RMS
         this->mRMS.setConstant(0.0);

         scalar_iterator_range sRange = this->scalarRange();
         assert(std::distance(sRange.first, sRange.second) == 1);

         // go through each vertical level and find the 0th mode (horizontal average) 
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
            // Mean
            this->mKurt(k_) = (this->rInVar.phys().slice(k) - mAvg(k_)).array().pow(4).sum()/(this->mRMS(k_)).pow(4);
         }


      }

      void Cartesian1DScalarKurtWriter::postcompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

#ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mKurt.data(), this->mKurt.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif //GEOMHDISCC_MPI
      }

      void Cartesian1DScalarKurtWriter::prewrite()
      {
         IVariableAsciiEWriter::prewrite();

         if(FrameworkMacro::allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mZ.transpose() << std::endl;
         }
      }

      void Cartesian1DScalarKurtWriter::write()
      {
         // Create file
         this->preWrite();


         // Check if the workflow allows IO to be performed
         if(FrameworkMacro::allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mKurt.transpose() << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic energy is NaN
         if(std::isnan(this->mKurt.sum()))
         {
            FrameworkMacro::abort(99);

            throw Exception("Kurtosis is NaN!");
         }
      }

   }
}
