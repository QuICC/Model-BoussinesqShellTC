/**
 * @file TransforStepsMacro.h
 * @brief Preprocessor macros to setup correct transform steps depending on geometry 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef TRANSFORMSTEPSMACRO_H
#define TRANSFORMSTEPSMACRO_H

#if defined GEOMHDISCC_SPATIALSCHEME_TTT || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_FFF 
   // include transform steps for 3D cartesian geometry
   #include "TransformConfigurators/Cartesian3DTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_AFT
   // include transform steps for the cylindrical annulus geometry
   #include "TransformConfigurators/AnnulusTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_CFT || defined GEOMHDISCC_SPATIALSCHEME_WFT
   // include transform steps for the whole cylinder geometry
   #include "TransformConfigurators/CylinderTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_SLFL || defined GEOMHDISCC_SPATIALSCHEME_SLFM
   // include transform steps for the spherical shell geometry
   #include "TransformConfigurators/ShellTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_BLFL || defined GEOMHDISCC_SPATIALSCHEME_BLFM || defined GEOMHDISCC_SPATIALSCHEME_WLFL || defined GEOMHDISCC_SPATIALSCHEME_WLFM
   // include transform steps for the whole sphere geometry
   #include "TransformConfigurators/SphereTransformSteps.hpp"

#elif defined GEOMHDISCC_SPATIALSCHEME_TT || defined GEOMHDISCC_SPATIALSCHEME_TF || defined GEOMHDISCC_SPATIALSCHEME_FF
   // include transform steps for 2D cartesian geometry
   #include "TransformConfigurators/Cartesian2DTransformSteps.hpp"

#endif // defined GEOMHDISCC_SPATIALSCHEME_TTT || defined GEOMHDISCC_SPATIALSCHEME_TFT || defined GEOMHDISCC_SPATIALSCHEME_TFF || defined GEOMHDISCC_SPATIALSCHEME_FFF

#endif // TRANSFORMSTEPSMACRO_H
