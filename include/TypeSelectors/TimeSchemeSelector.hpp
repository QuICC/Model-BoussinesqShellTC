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

// Configure code to use ImExRK3
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRK3

   #include "Timesteppers/ImExRK3.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRK3 IntegratorSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRK3

// Configure code to use ImExSBDF2
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

   #include "Timesteppers/ImExSBDF2.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExSBDF2 IntegratorSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

#endif // TIMESCHEMESELECTOR_HPP
