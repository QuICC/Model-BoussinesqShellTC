/** 
 * @file KineticEnergyBeta3DQGPerWriter.cpp
 * @brief Source of the implementation of the ASCII kinetic energy writer for the Beta 3DQG periodic model
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
#include "IoVariable/KineticEnergyBeta3DQGPerWriter.hpp"

// Project includes
//
#include "Enums/Dimensions.hpp"
#include "Enums/FieldIds.hpp"
#include "ScalarFields/FieldTools.hpp"
#include "IoVariable/KineticEnergyTags.hpp"
#include "IoTools/IdToHuman.hpp"
#include "Python/PythonWrapper.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   KineticEnergyBeta3DQGPerWriter::KineticEnergyBeta3DQGPerWriter(std::string type)
      : IVariableAsciiEWriter(KineticEnergyTags::BASENAME, KineticEnergyTags::EXTENSION, KineticEnergyTags::HEADER, type, KineticEnergyTags::VERSION, Dimensions::Space::SPECTRAL)
   {
   }

   KineticEnergyBeta3DQGPerWriter::~KineticEnergyBeta3DQGPerWriter()
   {
   }

   void KineticEnergyBeta3DQGPerWriter::write()
   {
      // Create file
      this->preWrite();

      KineticEnergyBeta3DQGPerWriter::scalar_iterator sIt;
      KineticEnergyBeta3DQGPerWriter::scalar_iterator_range sRange = this->scalarRange();
      assert(std::distance(sRange.first, sRange.second) == 3);

      Array energy = Array::Zero(3);
      ArrayI mode = sRange.first->second->dom(0).spRes()->cpu()->dim(Dimensions::Transform::TRA1D)->mode(0);
      if(mode(2) == 0 && mode(3) == 0)
      {
         for(sIt = sRange.first; sIt != sRange.second; ++sIt)
         {
            if(sIt->first == PhysicalNames::KINETIC_ENERGY)
            {   
               energy(0) = sIt->second->dom(0).perturbation().point(0,0,0).real();
            } else if(sIt->first == PhysicalNames::ZONAL_KINETIC_ENERGY)
            {
               energy(1) = sIt->second->dom(0).perturbation().point(0,0,0).real();
            } else if(sIt->first == PhysicalNames::NONZONAL_KINETIC_ENERGY)
            {
               energy(2) = sIt->second->dom(0).perturbation().point(0,0,0).real();
            } else
            {
               throw Exception("Kinetic energy calculation failed!");
            }
         }
      }

      // Get the "global" Kinetic energy from MPI code
      #ifdef GEOMHDISCC_MPI
         MPI_Allreduce(MPI_IN_PLACE, energy.data(), energy.size(), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      #endif //GEOMHDISCC_MPI

      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile << std::setprecision(14) << this->mTime << "\t" << energy(0) << "\t" << energy(1) << "\t" << energy(2) << std::endl;
      }

      // Close file
      this->postWrite();

      // Abort if kinetic energy is NaN
      if(std::isnan(energy))
      {
         #ifdef GEOMHDISCC_MPI
            MPI_Abort(MPI_COMM_WORLD, 99);
         #endif //GEOMHDISCC_MPI

         throw Exception("Kinetic energy is NaN!");
      }
   }

}
}
