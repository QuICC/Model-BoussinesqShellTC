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

// Configure code to use 2R implementation
#if defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB2 || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3A || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3B || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3C || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3D || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3E

   #include "Timesteppers/SparseImExRK2RTimestepper.hpp"

   // Workaround for template typedef
   #define TimeSchemeTypeSelector SparseImExRK2RTimestepper

#endif //defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB2 || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3A || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3B || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3C || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3D || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3E

// Configure code to use #R implementation
#if defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3F || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB4

   #include "Timesteppers/SparseImExRK3RTimestepper.hpp"

   // Workaround for template typedef
   #define TimeSchemeTypeSelector SparseImExRK3RTimestepper

#endif //defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB3F || defined GEOMHDISCC_TIMESTEPPER_IMEXRKCB4

// Configure code to use to use old SparseTimestepper
#if defined GEOMHDISCC_TIMESTEPPER_IMEXRK3 || defined GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

   #include "Timesteppers/SparseOldImExTimestepper.hpp"

   // Workaround for template typedef
   #define TimeSchemeTypeSelector SparseOldImExTimestepper

#endif //defined GEOMHDISCC_TIMESTEPPER_IMEXRK3 || defined GEOMHDISCC_TIMESTEPPER_IMEXSBDF2

#endif // TIMESCHEMETYPESELECTOR_HPP
