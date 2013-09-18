/** 
 * @file GxlWriter.hpp 
 * @brief Implementation of the GXL format file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef GXLWRITER_HPP
#define GXLWRITER_HPP

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
#include "IoXml/IGxlFile.hpp"

namespace GeoMHDiSCC {

namespace IoXml {

   /**
    * @brief Implementation of the GXL format file writer
    */
   class GxlWriter: public IGxlFile<IoXml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name Name of the file
          */
         GxlWriter(const std::string& name);

         /**
          * @brief Destructor
          */
         virtual ~GxlWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();

         void graphResolution(SharedResolution spRes);
         
      protected:
         void createAttr(rapidxml::xml_node<>* parent, const std::string& name, const std::string& value);
         void graph1DResolution(SharedResolution spRes);
         void graph2DResolution(SharedResolution spRes);
         void graph3DResolution(SharedResolution spRes);

      private:
   };

   /// Typedef for a smart pointer of a GxlWriter
   typedef SharedPtrMacro<GxlWriter> SharedGxlWriter;

}
}

#endif // GXLWRITER_HPP
