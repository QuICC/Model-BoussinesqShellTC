/** 
 * @file EquationEigenSelector.hpp
 * @brief Selector to define equation eigen tools
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef EQUATIONEIGENSELECTOR_HPP
#define EQUATIONEIGENSELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Equations/CouplingInformation.hpp"
#include "Equations/Tools/IEigenTools.hpp"

namespace GeoMHDiSCC {

namespace Equations {

   SharedIEigenTools eigenSelector(const CouplingInformation::IndexType indexType);

}
}

#endif // EQUATIONEIGENSELECTOR_HPP
