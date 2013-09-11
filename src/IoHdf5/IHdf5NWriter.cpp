/** 
 * @file IHdf5NWriter.cpp
 * @brief Source of the implementation of a numbering HDF5 file
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

// System includes
//
#include <sstream>
#include <iomanip>

// External includes
//

// Class include
//
#include "IoHdf5/IHdf5NWriter.hpp"

// Project includes
//

namespace GeoMHDiSCC {

namespace IoHdf5 {

   const int IHdf5NWriter::ZEROWIDTH = 4;

   IHdf5NWriter::IHdf5NWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : IHdf5Writer(name, ext, header, type, version), mCounter(0), mBaseName(name)
   {
   }

   IHdf5NWriter::~IHdf5NWriter()
   {
   }

   void IHdf5NWriter::changeBasename(std::string base)
   {
      this->mBaseName = base;
   }

   void IHdf5NWriter::updateName()
   {
      // Create stringstream
      std::ostringstream   oss;

      // Create zerofilled counter number
      oss << std::setfill('0') << std::setw(IHdf5NWriter::ZEROWIDTH) << this->mCounter;

      // Replace filenme with new name
      this->resetName(this->mBaseName + oss.str());
   }

   void IHdf5NWriter::init()
   {
   }

   void IHdf5NWriter::preWrite()
   {
      // Update filename
      this->updateName();
      
      // Create new file
      this->open();
   }

   void IHdf5NWriter::postWrite()
   {
      // Close the file
      this->close();

      // Increment counter
      ++this->mCounter;
   }

   void IHdf5NWriter::finalize()
   {
   }

}
}
