/**
 * @file IAsciiRWriter.hpp
 * @brief Interface to an overwriting ASCII writer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef IASCIIRWRITER_HPP
#define IASCIIRWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IoAscii/IAsciiWriter.hpp"

namespace GeoMHDiSCC {

namespace IoAscii {

   /**
    * @brief Interface to an overwriting ASCII writer
    */
   class IAsciiRWriter: public IAsciiWriter
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
         IAsciiRWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiRWriter();

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

      private:
   };
}
}

#endif // IASCIIRWRITER_HPP
