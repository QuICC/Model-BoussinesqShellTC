/** 
 * @file NusseltCubicWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer for a 3D finite box
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
#include "IoVariable/NusseltCubicWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   NusseltCubicWriter::NusseltCubicWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   NusseltCubicWriter::~NusseltCubicWriter()
   {
   }

   void NusseltCubicWriter::write()
   {
      // Create file
      this->preWrite();

      NusseltCubicWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      NusseltCubicWriter::scalar_iterator sit = sRange.first;


      // Get the "global" Nusselt number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << nusselt << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
