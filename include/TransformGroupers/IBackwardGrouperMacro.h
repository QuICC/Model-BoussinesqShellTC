/**
 * @file IBackwardGrouperMacro.h
 * @brief Preprocessor macros to setup correct backward grouper
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IBACKWARDGROUPERMACRO_H
#define IBACKWARDGROUPERMACRO_H

#ifdef QUICC_SPATIALDIMENSION_3D
   #include "TransformGroupers/IBackwardGrouper3D.hpp"
#else
   #include "TransformGroupers/IBackwardGrouper2D.hpp"
#endif // QUICC_SPATIALDIMENSION_3D

#endif // IBACKWARDGROUPERMACRO_H
