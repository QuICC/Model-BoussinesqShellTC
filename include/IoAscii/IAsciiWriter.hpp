/** \file IAsciiWriter.hpp
 *  \brief General interface to an ASCII writer
 *
 *  \mhdBug Needs test
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

namespace GeoMHDiSCC {

namespace IoAscii {

   /**
    * \brief General interface to an ASCII writer
    */
   class IAsciiWriter: public AsciiFile
   {
      public:
         /**
         * @brief Constructor
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of file
         * @param version  Version string of file 
         */
         IAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init() = 0;

         /**
          * @brief Write the content
          */
         virtual void write() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize() = 0;
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::ofstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

      private:
   };

   /// Typedef for a shared pointer of a IAsciiWriter
   typedef SharedPtrMacro<IAsciiWriter> SharedIAsciiWriter;
}
}

#endif // IASCIIWRITER_HPP
