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
         
         /**
          * @brief Create communication graph
          */
         void graphCommunication(const std::vector<std::multimap<int,int> >& structure);
         
      protected:
         /**
          * @brief Create an attr tag in the xml tree
          *
          * @param parent  Parent node to attach to
          * @param name    Name of the attr tag
          * @param value   Value to store in string child
          */
         void createAttr(rapidxml::xml_node<>* parent, const std::string& name, const std::string& value);

      private:
   };

   /// Typedef for a smart pointer of a GxlWriter
   typedef SharedPtrMacro<GxlWriter> SharedGxlWriter;

}
}

#endif // GXLWRITER_HPP
