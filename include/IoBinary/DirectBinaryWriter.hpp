/** 
 * @file DirectBinaryWriter.hpp
 * @brief Implementation of a direct access binary writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef DIRECTBINARYWRITER_HPP
#define DIRECTBINARYWRITER_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IoBinary/IBinaryWriter.hpp"

namespace QuICC {

namespace IoBinary {

   /**
    *@brief Implementation of a direct access binary writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
    */
   class DirectBinaryWriter: public IBinaryWriter
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
         DirectBinaryWriter(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~DirectBinaryWriter();

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

#endif // DIRECTBINARYWRITER_HPP
