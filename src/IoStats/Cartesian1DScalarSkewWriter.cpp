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
#include "IoStats/SkewTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

   namespace IoStats {

      Cartesian1DScalarSkewWriter::Cartesian1DScalarSkewWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type)
         : IStatisticsAsciiEWriter(prefix + SkewTags::BASENAME, SkewTags::EXTENSION, prefix + SkewTags::HEADER, type, SkewTags::VERSION, Dimensions::Space::SPECTRAL), mArea(-1), mAvg(Avg), mRMS(RMS), mSkew(-Array::Ones(2))
      {
      }

      Cartesian1DScalarSkewWriter::~Cartesian1DScalarSkewWriter()
      {
      }

      void Cartesian1DScalarSkewWriter::init()
      {

         IStatisticsAsciiEWriter::init();

         if(FrameworkMacro::allowsIO())
         {
            this->mFile << "# " << std::setprecision(14) << (1.0 + this->mMesh.at(0).transpose().reverse().array())/2.0 <<std::endl;
         }
         int dimZ = this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::PHYSICAL);
         int dimX = this->mspRes->sim()->dim(Dimensions::Simulation::SIM2D, Dimensions::Space::PHYSICAL);
         int dimY = this->mspRes->sim()->dim(Dimensions::Simulation::SIM3D, Dimensions::Space::PHYSICAL);

         this->mArea = dimX*dimY;
         this->mSkew = Array::Zero(dimZ);
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
            //this->mSkew(k_) = (sRange.first->second->dom(0).phys().slice(k).array() - mAvg->average()(k_)).array().pow(3).sum()/(mRMS->RMS()(k_)).pow(3);
            this->mSkew(k_) = (sRange.first->second->dom(0).phys().slice(k).array() - mAvg->average()(k_)).array().pow(3).sum();
            this->mSkew(k_) = this->mSkew(k_)/this->mArea;
         }

      }

      void Cartesian1DScalarSkewWriter::postCompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

#ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mSkew.data(), this->mSkew.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif //GEOMHDISCC_MPI

            this->mSkew = this->mSkew.array()/(mRMS->RMS().array()*mRMS->RMS().array()*mRMS->RMS().array());
      }


      void Cartesian1DScalarSkewWriter::write()
      {
         // Create file
         this->preWrite();


         // Check if the workflow allows IO to be performed
         if(FrameworkMacro::allowsIO())
         {
            this->mFile << std::setprecision(14) << this->mTime << " \t" << this->mSkew.transpose() << std::endl;
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
