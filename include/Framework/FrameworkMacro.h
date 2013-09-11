/**
 * @file FrameworkMacro.h
 * @brief Preprocessor macros to setup the framework depending on parallelisation options 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef FRAMEWORKMACRO_H
#define FRAMEWORKMACRO_H

#ifdef GEOMHDISCC_MPI
   // include MPI framework
   #include "Framework/MpiFramework.hpp"

   namespace GeoMHDiSCC {
      /// Typedef for a generic framework based on the MpiFramework
      typedef MpiFramework  FrameworkMacro;
   }
#else
   // include serial framework
   #include "Framework/SerialFramework.hpp"

   namespace GeoMHDiSCC {

      /// Typedef for a generic framework based on the SerialFramework
      typedef SerialFramework  FrameworkMacro;
   }
#endif // GEOMHDISCC_MPI

#endif // FRAMEWORKMACRO_H
