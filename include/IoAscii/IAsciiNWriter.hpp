/**
 * @file IAsciiNWriter.hpp
 * @brief Interface to a numbered file ASCII writer 
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef IASCIINWRITER_HPP
#define IASCIINWRITER_HPP

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
    * @brief Interface to a numbered file ASCII writer
    */
   class IAsciiNWriter: public IAsciiWriter
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
         IAsciiNWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IAsciiNWriter();

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
   };
}
}

#endif // IASCIINWRITER_HPP
