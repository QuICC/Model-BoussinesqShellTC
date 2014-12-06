/**
 * @file TransforStepsMacro.h
 * @brief Preprocessor macros to setup correct transform steps depending on geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMSTEPSMACRO_H
#define TRANSFORMSTEPSMACRO_H

#if defined GEOMHDISCC_SPATIALSCHEME_TTT || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_FFF 
   // include transform steps for cartesian geometry
   #include "TransformConfigurators/CartesianTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_CFT || defined GEOMHDISCC_SPATIALSCHEME_AFT || defined GEOMHDISCC_SPATIALSCHEME_WFT
   // include transform steps for cylindrical geometry
   #include "TransformConfigurators/CylindricalTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_SLF || defined GEOMHDISCC_SPATIALSCHEME_BLF || defined GEOMHDISCC_SPATIALSCHEME_WLF
   // include transform steps for spherical geometry
   #include "TransformConfigurators/SphericalTransformSteps.hpp"

#endif // defined GEOMHDISCC_SPATIALSCHEME_TTT || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_FFF

#endif // TRANSFORMSTEPSMACRO_H
