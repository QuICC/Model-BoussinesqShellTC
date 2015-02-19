/** 
 * @file VtpWriter.hpp 
 * @brief Implementation of the VTK PolyData XML format file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef VTPWRITER_HPP
#define VTPWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoXml/IXmlWriter.hpp"
#include "IoXml/IVtpFile.hpp"

namespace GeoMHDiSCC {

namespace IoXml {

   /**
    * @brief Implementation of the VTK PolyData XML format file writer
    */
   class VtpWriter: public IVtpFile<IoXml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         VtpWriter(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~VtpWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();

         /**
          * @brief Convert resolution to XML for VTP file
          */
         void representResolution(SharedTransformResolution spRes, const int rank);

      protected:

      private:
   };

   /// Typedef for a smart pointer of a VtpWriter
   typedef SharedPtrMacro<VtpWriter> SharedVtpWriter;

}
}

#endif // VTPWRITER_HPP
