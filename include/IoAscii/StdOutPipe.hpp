/**
 * @file StdOutPipe.hpp
 * @brief Implementation of standard output pipe file 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef STDOUTPIPE_HPP
#define STDOUTPIPE_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoAscii/IAsciiEWriter.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   /**
    * @brief Implementation of the standard output pipe into an ASCII file
    */
   class StdOutPipe: public IAsciiEWriter
   {
      public:
         /**
         * @brief Constructor 
         *
         * @param name Name of the std message ouput file
         */
         StdOutPipe(std::string name);

         /**
         * @brief Destructor
         */
         ~StdOutPipe();

         /**
          * @brief Init the file
          */
         void init();

         /**
          * @brief
          */
         virtual void write();

         /**
          * @brief Finalise
          */
         void finalize();
         
      protected:

      private:
         /**
          * @brief HEADER part for StdOutPipe file
          */
         static const std::string  HEADER;

         /**
          * @brief TYPE part for StdOutPipe file
          */
         static const std::string  TYPE;

         /**
          * @brief VERSION part for StdOutPipe file
          */
         static const std::string  VERSION;

         /**
          * @brief BASENAME of StdOutPipe file
          */
         static const std::string  BASENAME;

         /**
          * @brief EXTENSION of StdOutPipe file
          */
         static const std::string  EXTENSION;

         /**
          * Backup std::cout buffer
          */
         std::streambuf   *mpCoutBuffer;
   };

   /// Typedef for a smart shared pointer of a StdOutPipe
   typedef SharedPtrMacro<StdOutPipe> SharedStdOutPipe;

}
}

#endif // STDOUTPIPE_HPP
