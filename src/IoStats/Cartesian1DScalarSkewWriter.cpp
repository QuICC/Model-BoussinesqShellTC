/** 
 * @file Cartesian1DScalarSkewWriter.cpp
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
#include "IoStats/Cartesian1DScalarSkewWriter.hpp"

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

      Cartesian1DScalarSkewWriter::Cartesian1DScalarSkewWriter(const std::string& prefix, const std::string& type)
         : IVariableAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL), mEnergy(-Array::Ones(2))
      {
      }

      Cartesian1DScalarSkewWriter::~Cartesian1DScalarSkewWriter()
      {
      }

      void Cartesian1DScalarSkewWriter::init()
      {
         
         IVariableAsciiEWriter::init();
      }

      void Cartesian1DScalarSkewWriter::precompute(Transform::TransformCoordinatorType& coord)
      {
         // nothing here
      }

      void Cartesian1DScalarSkewWriter::compute(Transform::TransformCoordinatorType& coord)
      {
         // calculate the transforms and calculate the Skew

         scalar_iterator_range sRange = this->scalarRange();
         assert(std::distance(sRange.first, sRange.second) == 1);

         // go through each vertical level and find the difference from each point to the 0th mode (horizontal average) 
         for(int k = 0; k < this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->dim<Dimensions::Data::DAT3D>(); ++k)
         {
            int k_ = this->mspRes->cpu()->dim(Dimensions::Transform::TRA3D)->idx<Dimensions::Data::DAT3D>(k);
            // Skew calculation
            this->mSkew(k_) = (this->rInVar.phys().slice(k) - mAvg(k_)).array().pow(3).sum()/(this->mRMS(k_)).pow(3);
         }

      }

      void Cartesian1DScalarSkewWriter::postcompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

#ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mSkew.data(), this->mSkew.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif //GEOMHDISCC_MPI
      }

      void Cartesian1DScalarSkewWriter::prewrite()
      {
         IVariableAsciiEWriter::prewrite();

         if(FrameworkMacro::allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mZ.transpose() << std::endl;
         }
      }

      void Cartesian1DScalarSkewWriter::write()
      {
         // Create file
         this->preWrite();


         // Check if the workflow allows IO to be performed
         if(FrameworkMacro::allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mRMS.transpose() << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic energy is NaN
         if(std::isnan(this->mSkew.sum()))
         {
            FrameworkMacro::abort(99);

            throw Exception("Skew is NaN!");
         }
      }

   }
}
