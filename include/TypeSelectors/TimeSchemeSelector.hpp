/** 
 * @file TimeSchemeSelector.hpp
 * @brief Definition of some useful typedefs for the timestep scheme
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESCHEMESELECTOR_HPP
#define TIMESCHEMESELECTOR_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Timesteppers/ImExRK3.hpp"

namespace GeoMHDiSCC {

   namespace Timestep {

      typedef ImExRK3 TimeSchemeType;

   }
}

#endif // TIMESCHEMESELECTOR_HPP
