/** \file StateFileTags.cpp
 *  \brief Source of the definitions and names used by the state file readers/writers
 */

// System includes
//

// External includes
//

// Class include
//
#include "IoVariable/StateFileTags.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoVariable {

   const std::string StateFileTags::HEADER = "StateFile";

   const std::string StateFileTags::VERSION = "1.0";

   const std::string StateFileTags::BASENAME = "state";

   const std::string StateFileTags::EXTENSION = ".hdf5";
}
}