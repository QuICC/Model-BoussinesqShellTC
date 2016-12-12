/**
 * @file BackwardConfiguratorMacro.h
 * @brief Preprocessor macros to setup correct backward configurator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef BACKWARDCONFIGURATORMACRO_H
#define BACKWARDCONFIGURATORMACRO_H

#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformConfigurators/BackwardConfigurator3D.hpp"
#else
   #include "TransformConfigurators/BackwardConfigurator2D.hpp"
#endif // GEOMHDISCC_SPATIALDIMENSION_3D

#endif // BACKWARDCONFIGURATORMACRO_H
