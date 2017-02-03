/**
 * @file IAsciiReader.hpp
 * @brief Interface to an ASCII file reader 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IASCIIREADER_HPP
#define IASCIIREADER_HPP

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
    * @brief Interface to an ASCII file reader
    */
   class IAsciiReader: public AsciiFile
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
         IAsciiReader(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiReader();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Read the content
          */

         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::ifstream mFile;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Check compatibility of opened file
          */
         virtual void checkCompatibility();

      private:
   };

}
}

#endif // IASCIIREADER_HPP
