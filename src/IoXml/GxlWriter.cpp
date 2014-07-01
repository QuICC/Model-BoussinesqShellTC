/** 
 * @file GxlWriter.cpp
 * @brief Source of the implementation of the GXL format file writer
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//
#include <set>
#include <tr1/tuple>

// External includes
//
#include <rapidxml_print.hpp>

// Class include
//
#include "IoXml/GxlWriter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Enums/DimensionTools.hpp"

namespace GeoMHDiSCC {

namespace IoXml {

   GxlWriter::GxlWriter(const std::string& name)
      : IGxlFile<IoXml::IXmlWriter>(name)
   {
   }

   GxlWriter::~GxlWriter()
   {
   }

   void GxlWriter::graphCommunication(const std::vector<std::multimap<int,int> >& structure)
   {
      std::stringstream oss;

      // Set subgraph long identifier
      std::vector<std::string> lName;
      lName.push_back("cluster_A");
      lName.push_back("cluster_B");

      // Set subgraph short identifier
      std::vector<std::string> sName;
      sName.push_back("A");
      sName.push_back("B");

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the two transposes
      for(size_t i = 0; i < structure.size(); i++)
      {
         // Create subgraph
         rapidxml::xml_node<> *pNSubgraph = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "N_" << lName.at(i);
         pNSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
         oss.str("");
         rapidxml::xml_node<> *pSubgraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
         pSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(lName.at(i).c_str(),0)));

         if(static_cast<Dimensions::Transform::Id>(i) == Dimensions::Transform::TRA1D)
         {
            this->createAttr(pSubgraph, "label", "Transpose 1D/2D");
         } else
         {
            this->createAttr(pSubgraph, "label", "Transpose 2D/3D");
         }
         this->createAttr(pSubgraph, "color", color.at(i));

         pNSubgraph->append_node(pSubgraph);
         pGraph->append_node(pNSubgraph);

         // Loop over the CPUs
         for(std::multimap<int,int>::const_iterator itCpu = structure.at(i).begin(); itCpu != structure.at(i).end(); itCpu = structure.at(i).upper_bound(itCpu->first))
         {
            rapidxml::xml_node<> *pCpu = this->mXML.allocate_node(rapidxml::node_element, "node");
            oss << "cpu" << itCpu->first << sName.at(i);
            pCpu->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            oss << "CPU " << itCpu->first;
            this->createAttr(pCpu, "label", oss.str());
            oss.str("");

            pSubgraph->append_node(pCpu);
         }

         // Loop over the CPUs
         for(std::multimap<int,int>::const_iterator itCpu = structure.at(i).begin(); itCpu != structure.at(i).end(); ++itCpu)
         {
            // Set "from" attribute
            rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
            oss << "cpu" << itCpu->first << sName.at(i);
            pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            // Set "to" attribute
            oss << "cpu" << itCpu->second << sName.at(i);
            pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            // Add edge
            this->createAttr(pEdge, "color", color.at(i));
            pSubgraph->append_node(pEdge);
         }
      }
   }

   void GxlWriter::createAttr(rapidxml::xml_node<>* sParent, const std::string& name, const std::string& value)
   {
      rapidxml::xml_node<>* pAttr;
      rapidxml::xml_node<>* pStr;

      pAttr = this->mXML.allocate_node(rapidxml::node_element, "attr");
      pAttr->append_attribute(this->mXML.allocate_attribute("name", this->mXML.allocate_string(name.c_str()),0));
      pStr = this->mXML.allocate_node(rapidxml::node_element, "string");
      pStr->value(this->mXML.allocate_string(value.c_str(),0));

      pAttr->append_node(pStr);
      sParent->append_node(pAttr);
   }

   void GxlWriter::write()
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
}
}
