/** 
 * @file VisualizationFileTags.cpp
 * @brief Source of the definitions and names used by the visualisation file writers
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/VisualizationFileTags.hpp"

// Project includes
//

namespace QuICC {

namespace IoVariable {

   const std::string VisualizationFileTags::HEADER = "VisualizationFile";

   const std::string VisualizationFileTags::VERSION = "1.0";

   const std::string VisualizationFileTags::BASENAME = "visState";

   const std::string VisualizationFileTags::EXTENSION = ".hdf5";

   const std::string VisualizationFileTags::MESH = "mesh";

   const std::string VisualizationFileTags::GRID = "grid";
}
}
