/** 
 * @file Cartesian1DNusseltZWriter.cpp
 * @brief Source of the implementation of the ASCII Nusselt number writer through the Z boundary
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
#include "IoVariable/Cartesian1DNusseltZWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/NusseltTags.hpp"
#include "IoTools/IdToHuman.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   Cartesian1DNusseltZWriter::Cartesian1DNusseltZWriter(std::string type)
      : IVariableAsciiEWriter(NusseltTags::BASENAME, NusseltTags::EXTENSION, NusseltTags::HEADER, type, NusseltTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   Cartesian1DNusseltZWriter::~Cartesian1DNusseltZWriter()
   {
   }

   void Cartesian1DNusseltZWriter::write()
   {
      // Create file
      this->preWrite();

      Cartesian1DNusseltZWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
      Cartesian1DNusseltZWriter::scalar_iterator sit = sRange.first;

      ArrayI mode = sit->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      MHDFloat nusselt = 0.0;
      if(mode(2) == 0 && mode(3) == 0)
      {
         // Create boundary operator
         Array bc = 2.0*Array::Ones(this->mspRes->sim()->dim(Dimensions::Simulation::SIM1D, Dimensions::Space::SPECTRAL));
         bc(0) = bc(0)/2.0;

         // Compute Nusselt number
         nusselt = -bc.dot(sit->second->dom(0).perturbation().profile(0,0).real());
      }

      // Get the "global" Nusselt number from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &nusselt, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << nusselt << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if Nusselt number is NaN
      if(std::isnan(nusselt))
      {
         FrameworkMacro::abort(99);

         throw Exception("Nusselt number is NaN!");
      }
   }

}
}
