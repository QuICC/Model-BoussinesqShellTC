/** \file IAsciiEWriter.hpp
 *  \brief General interface of an "extending" ASCII writer
 */

#ifndef IASCIIEWRITER_HPP
#define IASCIIEWRITER_HPP

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
    * \brief General interface of an "extending" ASCII writer
    */
   class IAsciiEWriter: public IAsciiWriter
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
         IAsciiEWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiEWriter() {};

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
         virtual void finalise();
         
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

#endif // IASCIIEWRITER_HPP