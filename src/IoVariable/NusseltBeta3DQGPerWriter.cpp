/** 
 * @file NusseltBeta3DQGPerWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer for the Beta 3DQG periodic model
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
#include "IoVariable/NusseltBeta3DQGPerWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   NusseltBeta3DQGPerWriter::NusseltBeta3DQGPerWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   NusseltBeta3DQGPerWriter::~NusseltBeta3DQGPerWriter()
   {
   }

   void NusseltBeta3DQGPerWriter::write()
   {
      // Create file
      this->preWrite();

      NusseltBeta3DQGPerWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      NusseltBeta3DQGPerWriter::scalar_iterator sit = sRange.first;

      MHDComplex nusselt = MHDComplex(0.0,0.0);
      for(int i = 0; i < sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->dim<Dimensions::Data::DAT3D>(); ++i)
      {
         if(sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT2D>(0,i) == 0)
         {
            MHDComplex zi_ = MHDComplex(0.0, static_cast<MHDFloat>(sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->idx<Dimensions::Data::DAT3D>(i)));
            nusselt += zi_*sit->second->dom(0).perturbation().point(0,0,i);
         }
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << 1.0 + nusselt.real() << "\t" << nusselt.imag() << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
