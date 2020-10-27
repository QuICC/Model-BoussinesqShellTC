/**
 * @file IAsciiWriter.hpp
 * @brief General interface to an ASCII writer 
 */

#ifndef IASCIIWRITER_HPP
#define IASCIIWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <fstream>

// External includes
//

// Project includes
//
#include "IoAscii/AsciiFile.hpp"

namespace QuICC {

namespace IoAscii {

   /**
    * @brief General interface to an ASCII writer
    */
   class IAsciiWriter: public AsciiFile
   {
      public:
         /**
          * @brief Possible write modes
          */
         enum WriteMode {
            EXTEND, // New values are appended
            OVERWRITE, // File is overwritten
            NUMBER, // New file is created with incremented number
         };

         /**
         * @brief Constructor
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of file
         * @param version  Version string of file 
         * @param mode     Write mode of file
         */
         IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const WriteMode mode = EXTEND);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();
         
      protected:
         /**
          * @brief Operation to perform just before writing data
          */
         void preWrite();

         /**
          * @brief Operation to perform just after writing data
          */
         void postWrite();

         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Open the file in debug mode (no IO filter)
          */
         void openDebug();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Close the file in debug mode (no IO filter)
          */
         void closeDebug();

         /**
          * @brief Create the file
          */
         void create();

         /**
          * @brief end the file
          */
         void end();

         /**
          * @brief Update the file name with counter value
          */
         void updateName();

      private:
         /**
          * @brief Zero Fill width
          */
         static const int msIDWidth;

         /**
          * @brief File counter
          */
         int mCounter;

         /**
          * @brief Base name of the file
          */
         const std::string    mBaseName;

         /**
          * @brief Write mode of file
          */
         const WriteMode mMode;
   };

   /// Typedef for a shared pointer of a IAsciiWriter
   typedef SharedPtrMacro<IAsciiWriter> SharedIAsciiWriter;
}
}

#endif // IASCIIWRITER_HPP
