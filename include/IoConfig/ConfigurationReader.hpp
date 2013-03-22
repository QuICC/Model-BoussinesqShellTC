/** \file ConfigurationReader.hpp 
 *  \brief Implementation of the XML configuration file reader
 */

#ifndef CONFIGURATIONREADER_HPP
#define CONFIGURATIONREADER_HPP

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
#include "IoXml/IXmlReader.hpp"
#include "IoConfig/IConfigurationFile.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Implementation of the XML configuration file reader
    */
   class ConfigurationReader: public IConfigurationFile<IoXml::IXmlReader>
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param type Type of the simulation
          */
         ConfigurationReader(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ConfigurationReader();

         /**
          * @brief Read content of configuration file
          */
         virtual void read();
         
      protected:

      private:
         /**
          * @brief NAME of the configuration file
          */
         static const std::string NAME;
   };

   /// Typedef for a smart pointer of a ConfigurationReader
   typedef SharedPtrMacro<ConfigurationReader> SharedConfigurationReader;

}
}

#endif // CONFIGURATIONREADER_HPP
