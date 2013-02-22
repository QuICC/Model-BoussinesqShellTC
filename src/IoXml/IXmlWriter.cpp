/** \file IXmlWriter.cpp 
 *  \brief Source of the implementation of the XML writer
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//
#include <rapidxml_print.hpp>

// Class include
//
#include "IoXml/IXmlWriter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoXml {

   IXmlWriter::IXmlWriter(std::string name, std::string ext, std::string header, std::string type, std::string version)
      : XmlFile(name, ext, header, type, version)
   {
   }

   void IXmlWriter::init()
   {
      // Node variables
      rapidxml::xml_node<>* node;
      rapidxml::xml_node<>* child;

      // Create xml declaration
      node = this->mXML.allocate_node(rapidxml::node_declaration);
      node->append_attribute(this->mXML.allocate_attribute("version", "1.0"));
      node->append_attribute(this->mXML.allocate_attribute("encoding", "utf-8"));
      this->mXML.append_node(node);
      
      // FILEMETA node
      node = this->mXML.allocate_node(rapidxml::node_element, this->fileTag().c_str());
      this->mXML.append_node(node);
      
      // HEADER node
      child = this->mXML.allocate_node(rapidxml::node_element, this->headerTag().c_str());
      child->value(this->header().c_str());
      node->append_node(child);
      
      // TYPE node
      child = this->mXML.allocate_node(rapidxml::node_element, this->typeTag().c_str());
      child->value(this->type().c_str());
      node->append_node(child);
      
      // VERSION node
      child = this->mXML.allocate_node(rapidxml::node_element, this->versionTag().c_str());
      child->value(this->version().c_str());
      node->append_node(child);
   }

   void IXmlWriter::open()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Get handle to file
         this->mFile.open(this->filename().c_str());

         // Check that opening was a success
         if(! this->mFile.is_open())
         {
            throw Exception("Couldn't open XML file " + this->filename() + "!");
         }
      }
   }

   void IXmlWriter::finalise()
   {
   }

   void IXmlWriter::close()
   {
      // Check if the workflow allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         this->mFile.close();
      }
   }

   void IXmlWriter::write()
   {
      // Do pre writing processing
      this->preWrite();

      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Write xml content to file
         this->mFile << this->mXML;
      }

      // Do post writing processing
      this->postWrite();
   }

   void IXmlWriter::preWrite()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create the file
         this->open();
      }
   }

   void IXmlWriter::postWrite()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Close the file
         this->close();
      }
   }

}
}
