/** 
 * @file SphericalTorPolEnergyWriter.cpp
 * @brief Source of the implementation of the ASCII spherical harmonics energy calculation for scalar field
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
#include "IoVariable/SphericalTorPolEnergyWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "IoVariable/EnergyTags.hpp"
#include "TypeSelectors/TransformSelector.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   SphericalTorPolEnergyWriter::SphericalTorPolEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableHeavyAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   SphericalTorPolEnergyWriter::~SphericalTorPolEnergyWriter()
   {
   }

   void SphericalTorPolEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
   }

   void SphericalTorPolEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // get iterator to field
      vector_iterator_range vRange = this->vectorRange();
      assert(std::distance(vRange.first, vRange.second) == 1);
      assert(FieldComponents::Spectral::ONE == FieldComponents::Spectral::TOR);
      assert(FieldComponents::Spectral::TWO == FieldComponents::Spectral::POL);

      MHDFloat energy = 0.0;

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, &energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(16) << this->mTime << "\t" << energy << std::endl;
      }

      // Close file
      this->postWrite();
   }

}
}
