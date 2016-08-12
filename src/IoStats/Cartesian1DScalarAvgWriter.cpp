/** 
 * @file Cartesian1DScalarAvgWriter.cpp
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
#include "IoStats/Cartesian1DScalarAvgWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoStats/AvgTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

   namespace IoStats {

      Cartesian1DScalarAvgWriter::Cartesian1DScalarAvgWriter(const std::string& prefix, const std::string& type)
         : IStatisticsAsciiEWriter(prefix + AvgTags::BASENAME, AvgTags::EXTENSION, prefix + AvgTags::HEADER, type, AvgTags::VERSION, Dimensions::Space::SPECTRAL), mArea(-1), mAvg(-Array::Ones(2))
      {
      }

      Cartesian1DScalarAvgWriter::~Cartesian1DScalarAvgWriter()
      {
      }

      void Cartesian1DScalarAvgWriter::init()
      {
         IStatisticsAsciiEWriter::init();

         if(FrameworkMacro::allowsIO())
         {
            this->mFile << "# " << std::setprecision(14) <<  (1.0 + this->mMesh.at(0).transpose().array())/2.0 << std::endl;
         }
      }

      void Cartesian1DScalarAvgWriter::preCompute(Transform::TransformCoordinatorType& coord)
      {
         // Dealias variable data
         scalar_iterator_range sRange = this->scalarRange();

         if(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(0) == 0 && this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,0) == 0)
         {
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
         } else
         {
            this->mAvg = Array::Zero(this->mspRes->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DATF1D>());
         }
      }

      const Array& Cartesian1DScalarAvgWriter::average() const
      {
         return this->mAvg;
      }

      void Cartesian1DScalarAvgWriter::compute(Transform::TransformCoordinatorType& coord)
      {
         // nothing here for this one
      }

      void Cartesian1DScalarAvgWriter::postCompute(Transform::TransformCoordinatorType& coord)
      {
         // MPI gathering

         #ifdef GEOMHDISCC_MPI
            MPI_Allreduce(MPI_IN_PLACE, this->mAvg.data(), this->mAvg.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
         #endif //GEOMHDISCC_MPI
      }

      void Cartesian1DScalarAvgWriter::write()
      {
         // Create file
         this->preWrite();

         // Check if the workflow allows IO to be performed
         if(FrameworkMacro::allowsIO())
         {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << this->mAvg.transpose() << std::endl;
         }

         // Close file
         this->postWrite();

         // Abort if kinetic energy is NaN
         if(std::isnan(this->mAvg.sum()))
         {
            FrameworkMacro::abort(99);

            throw Exception("Horizontal Avg is NaN!");
         }
      }

   }
}
