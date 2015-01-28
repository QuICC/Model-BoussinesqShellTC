/** 
 * @file SphericalScalarEnergyWriter.cpp
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
#include "IoVariable/SphericalScalarEnergyWriter.hpp"

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

   SphericalScalarEnergyWriter::SphericalScalarEnergyWriter(const std::string& prefix, const std::string& type)
      : IVariableHeavyAsciiEWriter(prefix + EnergyTags::BASENAME, EnergyTags::EXTENSION, prefix + EnergyTags::HEADER, type, EnergyTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   SphericalScalarEnergyWriter::~SphericalScalarEnergyWriter()
   {
   }

   void SphericalScalarEnergyWriter::init()
   {
   }

   void SphericalScalarEnergyWriter::compute(Transform::TransformCoordinatorType& coord)
   {
      scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 1);
//
//      // Copy data to transform storage
//      bwdData.topRows(spSetup->specSize()) = sRange.first->dom(0).perturbation().comp(FieldComponents::Spectral::SCALAR).data();
//
//      // Compute backward transform
//      transform.project(fwdData, bwdData, TransformType::ProjectorType::PROJ, Arithmetics::SET);
//
//      // Compute quadratic term
//      fwdData = fwdData.array()*fwdData.conjugate().array();
//
//      // Compute forward transform
//      transform.integrate(bwdData, fwdData, TransformType::IntegratorType::INTG, Arithmetics::SET);
   }

   void SphericalScalarEnergyWriter::write()
   {
      // Create file
      this->preWrite();

      // Compute Chebyshev integral
      MHDFloat energy = 0;

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
