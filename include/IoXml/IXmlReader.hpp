/** 
 * @file IXmlReader.hpp
 * @brief Interface to an XML reader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 * @version 0.9.0
 * @date 2013-09-11
 */

#ifndef IXMLREADER_HPP
#define IXMLREADER_HPP

// System includes
//
#include <fstream>
#include <sstream>
#include <vector>

// External includes
//
#include "rapidxml.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "IoXml/XmlFile.hpp"

namespace GeoMHDiSCC {

namespace IoXml {

   /**
    * @brief Interface to and XML reader (based on rapidXML)
    */
   class IXmlReader: public XmlFile
   {
      public:
         /**
         * @brief Constructor 
         *
         * @param name     Filename
         * @param ext      File extension
         * @param header   Header string of file
         * @param type     Type string of the file
         * @param version  Version string of the file
         */
         IXmlReader(std::string name, std::string ext, std::string header, std::string type, std::string version);

         /**
         * @brief Destructor
         */
         virtual ~IXmlReader();

         /**
          * @brief Initialise the file
          */
         virtual void init();

         /**
          * @brief Read the content
          */
         virtual void read() = 0;

         /**
          * @brief Finalise the file
          */
         virtual void finalize();
         
      protected:
         /**
          * @brief Handle to the file
          */
         std::ifstream mFile;

         /**
          * @brief Content of file converted into a vector
          */
         std::vector<char> mContent;

         /**
          * @brief Open the file
          */
         void open();

         /**
          * @brief Close the file
          */
         void close();

         /**
          * @brief Parse the XML content
          */
         void parse();

         /**
          * @brief Check compatibility of opened file
          */
         void checkCompatibility();

         /**
          * @brief templated function to read value from XML tree
          *
          * rapidXML reads the value from the file as a string. The conversion to an other format is selected
          * by the template parameter (ie. type of val). If wrong type is given the obtained value will be wrong.
          *
          * @param val  Storage for read value
          * @param base XML tree node
          * @param tag  XML tag to get value
          *
          * \tparam T Type of the value to read
          */
         template <typename T> void readValue(T& val, rapidxml::xml_node<>* base, const std::string& tag);

      private:
   };

   template <typename T>  void IXmlReader::readValue(T& val, rapidxml::xml_node<>* base, const std::string& tag)
   {
      // Create string stream to do type conversion
      std::istringstream   iss;
      std::string tmp;

      // Extract node for given tag
      rapidxml::xml_node<> *vnode = base->first_node(tag.c_str());

      // Check for successful node extraction
      if(vnode)
      {
         // Get node value
         tmp = vnode->value();

         // Convert node value string into template type
         iss.str(tmp);
         iss >> val;
         iss.clear();

      } else
      {
         throw Exception("Unknown XML tag " + tag + "!");
      }
   }

}
}

#endif // IXMLREADER_HPP
