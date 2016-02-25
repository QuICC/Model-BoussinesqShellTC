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

// Configure code to use ImExRKCB2
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB2

   #include "Timesteppers/ImExRKCB2.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB2 TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB2

// Configure code to use ImExRKCB3a
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB3A

   #include "Timesteppers/ImExRKCB3a.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB3a TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB3A

// Configure code to use ImExRKCB3b
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB3B

   #include "Timesteppers/ImExRKCB3b.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB3b TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB3B

// Configure code to use ImExRKCB3c
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB3C

   #include "Timesteppers/ImExRKCB3c.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB3c TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB3C

// Configure code to use ImExRKCB3d
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB3D

   #include "Timesteppers/ImExRKCB3d.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB3d TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB3D

// Configure code to use ImExRKCB3e
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRKCB3E

   #include "Timesteppers/ImExRKCB3e.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRKCB3e TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB3E

// Configure code to use ImExRK3
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXRK3

   #include "Timesteppers/ImExRK3.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExRK3 TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXRK3

// Configure code to use ImExSBDF2
#ifdef GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

   #include "Timesteppers/ImExSBDF2.hpp"

   namespace GeoMHDiSCC {

      namespace Timestep {

         typedef ImExSBDF2 TimeSchemeSelector;

      }
   }
#endif //GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

#endif // TIMESCHEMESELECTOR_HPP
