/**
 * @file IForwardGrouperMacro.h
 * @brief Preprocessor macros to setup correct forward grouper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IFORWARDGROUPERMACRO_H
#define IFORWARDGROUPERMACRO_H

#ifdef GEOMHDISCC_SPATIALDIMENSION_3D
   #include "TransformGroupers/IForwardGrouper3D.hpp"
#else
   #include "TransformGroupers/IForwardGrouper2D.hpp"
#endif // GEOMHDISCC_SPATIALDIMENSION_3D

#endif // IFORWARDGROUPERMACRO_H