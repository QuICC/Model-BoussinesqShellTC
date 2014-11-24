/**
 * @file TransforStepsMacro.h
 * @brief Preprocessor macros to setup correct transform steps depending on geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMSTEPSMACRO_H
#define TRANSFORMSTEPSMACRO_H

#if defined GEOMHDISCC_CARTESIAN
   // include transform steps
   #include "TransformConfigurators/CartesianTransformSteps.hpp"

#elif defined GEOMHDISCC_CYLINDRICAL
   // include transform steps
   #include "TransformConfigurators/CartesianTransformSteps.hpp"

#elif defined GEOMHDISCC_SPHERICAL
   // include transform steps
   #include "TransformConfigurators/CartesianTransformSteps.hpp"

#endif // defined GEOMHDISCC_CARTESIAN

#endif // TRANSFORMSTEPSMACRO_H
