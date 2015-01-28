/** 
 * @file SimulationIoTools.cpp
 * @brief Source of the implementation of the tools for IO related calculations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "Simulation/SimulationIoTools.hpp"

// Project includes
//

namespace GeoMHDiSCC {

   SimulationIoTools::SimulationIoTools()
   {
   }

   SimulationIoTools::~SimulationIoTools()
   {
   }

   void SimulationIoTools::updateHeavyAscii(SimulationIoTools::ascii_iterator asciiBegin,  SimulationIoTools::ascii_iterator asciiEnd, Transform::TransformCoordinatorType& coord)
   {
      ascii_iterator it;
      for(it = asciiBegin; it != asciiEnd; ++it)
      {
         updateHeavyFile(*it, coord);
      }
   }

   void SimulationIoTools::updateHeavyFile(IoVariable::SharedIVariableHeavyAsciiEWriter spAscii, Transform::TransformCoordinatorType& coord)
   {
      spAscii->compute(coord);
   }
}
