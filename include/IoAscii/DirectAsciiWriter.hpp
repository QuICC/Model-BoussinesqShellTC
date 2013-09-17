/**
 * @file DirectAsciiWriter.hpp
 * @brief Implementation of a direct access ASCII writer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DIRECTASCIIWRITER_HPP
#define DIRECTASCIIWRITER_HPP

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
    * @brief Implementation of a direct access ASCII writer
    */
   class DirectAsciiWriter: public IAsciiWriter
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
         DirectAsciiWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~DirectAsciiWriter();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief This call does nothing in this case
          */
         virtual void write() {};

         /**
          * @brief Finalise the file
          */
         virtual void finalize();

         /**
          * @brief Get direct access to file handle
          */
         std::ofstream& file();
         
      protected:

      private:
   };
}
}

#endif // DIRECTASCIIWRITER_HPP
