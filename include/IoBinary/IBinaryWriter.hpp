/** \file IBinaryWriter.hpp
 *  \brief Interface to a general binary writer
 *
 *  \mhdBug Needs test
 */

#ifndef IBINARYWRITER_HPP
#define IBINARYWRITER_HPP

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
#include "IoBinary/BinaryFile.hpp"

namespace GeoMHDiSCC {

namespace IoBinary {

   /**
    * \brief Interface to a general binary writer
    */
   class IBinaryWriter: public BinaryFile
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
         IBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IBinaryWriter();

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

   /// Typedef for a shared pointer of a IBinaryWriter
   typedef SharedPtrMacro<IBinaryWriter> SharedIBinaryWriter;
}
}

#endif // IBINARYWRITER_HPP
