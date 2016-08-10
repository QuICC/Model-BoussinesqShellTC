/** 
 * @file Cartesian1DScalarRMSWriter.cpp
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
#include "IoStats/Cartesian1DScalarRMSWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoStats/RmsTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

   namespace IoStats {

      Cartesian1DScalarRMSWriter::Cartesian1DScalarRMSWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const std::string& type)
         : IStatisticsAsciiEWriter(prefix + RmsTags::BASENAME, RmsTags::EXTENSION, prefix + RmsTags::HEADER, type, RmsTags::VERSION, Dimensions::Space::SPECTRAL), mArea(-1), mAvg(Avg),mRMS(-Array::Ones(2))
      {
      }

      Cartesian1DScalarRMSWriter::~Cartesian1DScalarRMSWriter()
      {
      }

      void Cartesian1DScalarRMSWriter::init()
      {

         IStatisticsAsciiEWriter::init();

         if(FrameworkMacro::allowsIO())
         {
            this->mFile << "# " << std::setprecision(14) << this->mMesh.at(0).transpose() <<std::endl;
         }
         // Here is where we get the mean from the 0th mode of the variable in spectral space 

         int cols = this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>();

         // Normalize by Cartesian Area A = Nx*Ny (need to worry about 2 pi in fft?)
         this->mArea = cols*cols;

      }
      
      //RMS() function returns the array of RMS(Z) values when called by the Kurt and Skew functions
      const Array& Cartesian1DScalarRMSWriter::RMS() const
      {
         return this->mRMS;
      }

      void Cartesian1DScalarRMSWriter::compute(Transform::TransformCoordinatorType& coord)
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
            this->mRMS(k_) = (this->rInVar.phys().slice(k) - mAvg->average()(k_)).array().pow(2).sum();
         }


         // Normalize by the Cartesian volume
         this->mRMS /= this->mArea;

         //take sqroot
         this->mRMS = (this->mRMS).sqrt();
      }

      void Cartesian1DScalarRMSWriter::postCompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

#ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, this->mRMS.data(), this->mRMS.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif //GEOMHDISCC_MPI
      }


      void Cartesian1DScalarRMSWriter::write()
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
         if(std::isnan(this->mRMS.sum()))
         {
            FrameworkMacro::abort(99);

            throw Exception("RMS is NaN!");
         }
      }

   }
}
