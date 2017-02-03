/** 
 * @file VariableHdf5TagsDefs.cpp
 * @brief Source of the definitions and names used by the variable data files
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/VariableHdf5Tags.hpp"

// Project includes
//

namespace QuICC {

namespace IoVariable {

   const std::string VariableHdf5Tags::TRUNCATION = "truncation";

   const std::string VariableHdf5Tags::TRUNCPHYSICAL = "physical";

   const std::string VariableHdf5Tags::TRUNCTRANSFORM = "transform";

   const std::string VariableHdf5Tags::TRUNCSPECTRAL = "spectral";

   const std::string VariableHdf5Tags::TRUNCDIM = "dim";

   const std::string VariableHdf5Tags::PHYSICAL = "physical";

   const std::string VariableHdf5Tags::RUN = "run";

   const std::string VariableHdf5Tags::RUNTIME = "time";

   const std::string VariableHdf5Tags::RUNSTEP = "timestep";
}
}
