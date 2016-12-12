/**
 * @file ForwardConfiguratorMacro.h
 * @brief Preprocessor macros to setup correct forward configurator
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FORWARDCONFIGURATORMACRO_H
#define FORWARDCONFIGURATORMACRO_H

#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformConfigurators/ForwardConfigurator3D.hpp"
#else
   #include "TransformConfigurators/ForwardConfigurator2D.hpp"
#endif // GEOMHDISCC_SPATIALDIMENSION_3D

#endif // FORWARDCONFIGURATORMACRO_H
