/** 
 * @file TimeSchemeTypeSelector.hpp
 * @brief Definition of some useful typedefs for the timestep scheme type
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TIMESCHEMETYPESELECTOR_HPP
#define TIMESCHEMETYPESELECTOR_HPP

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

   #include "Timesteppers/SparseImExRK2RTimestepper.hpp"

   // Workaround for template typedef
   #define TimeSchemeTypeSelector SparseImExRK2RTimestepper

#endif //GEOMHDISCC_TIMESTEPPER_IMEXRKCB2

// Configure code to use to use old SparseTimestepper
#if defined GEOMHDISCC_TIMESTEPPER_IMEXRK3 || defined GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

   #include "Timesteppers/SparseTimestepper.hpp"

   // Workaround for template typedef
   #define TimeSchemeTypeSelector SparseTimestepper

#endif //defined GEOMHDISCC_TIMESTEPPER_IMEXRK3 || defined GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

#endif // TIMESCHEMETYPESELECTOR_HPP
