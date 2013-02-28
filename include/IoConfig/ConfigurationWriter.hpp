/** \file ConfigurationWriter.hpp 
 *  \brief Implementation of the XML configuration file writer
 *
 *  \mhdBug Needs test
 */

#ifndef CONFIGURATIONWRITER_HPP
#define CONFIGURATIONWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "IoXml/IXmlWriter.hpp"
#include "IoConfig/IConfigurationFile.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Implementation of the XML configuration file writer
    */
   class ConfigurationWriter: public IConfigurationFile<IoXml::IXmlWriter>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param type Type of the simulation
          */
         ConfigurationWriter(const int dim, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationWriter();

         /**
          * @brief Read content of configuration file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief NAME of the configuration file
          */
         static const std::string NAME;
   };

   /// Typedef for a smart pointer of a ConfigurationWriter
   typedef SharedPtrMacro<ConfigurationWriter> SharedConfigurationWriter;

}
}

#endif // CONFIGURATIONWRITER_HPP
