/**
 * @file FrameworkMacro.h
 * @brief Preprocessor macros to setup the framework depending on parallelisation options 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef FRAMEWORKMACRO_H
#define FRAMEWORKMACRO_H

#ifdef QUICC_MPI
   // include MPI framework
   #include "Framework/MpiFramework.hpp"

   namespace QuICC {
      /// Typedef for a generic framework based on the MpiFramework
      typedef MpiFramework  FrameworkMacro;
   }
#else
   // include serial framework
   #include "Framework/SerialFramework.hpp"

   namespace QuICC {

      /// Typedef for a generic framework based on the SerialFramework
      typedef SerialFramework  FrameworkMacro;
   }
#endif // QUICC_MPI

#endif // FRAMEWORKMACRO_H
