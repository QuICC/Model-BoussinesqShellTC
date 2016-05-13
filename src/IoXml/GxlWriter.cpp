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
#include "IoTools/IdToHuman.hpp"

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

   void GxlWriter::graphTransformPath(const std::map<PhysicalNames::Id, std::vector<Transform::TransformPath> >& paths, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");
      color.push_back("red");

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the paths
      int unid = 0;
      for(std::map<PhysicalNames::Id,std::vector<Transform::TransformPath> >::const_iterator nameIt = paths.begin(); nameIt != paths.end(); ++nameIt)
      {
         for(std::vector<Transform::TransformPath>::const_iterator pathIt = nameIt->second.begin(); pathIt != nameIt->second.end(); ++pathIt)
         {
            std::string field;

            rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
            oss << "p" << unid << "_" << nameIt->first << "_" << pathIt->startId();
            field = oss.str();
            oss.str("");
            pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(field.c_str(),0)));

            oss << IoTools::IdToHuman::toString(nameIt->first) << " ";
            if(dir == TransformDirection::FORWARD)
            {
               oss << IoTools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(pathIt->startId()));
            } else
            {
               oss << IoTools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(pathIt->startId()));
            }
            this->createAttr(pNode, "label", oss.str());
            oss.str("");

            pGraph->append_node(pNode);

            // Add all path nodes
            std::string curNode = field;
            for(int i = 0; i < pathIt->nEdges(); ++i)
            {
               rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
               oss << field << "_" << i;
               pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));

               // Set "from" attribute
               rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
               pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(curNode.c_str(),0)));

               // Set "to" attribute
               pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(oss.str().c_str(),0)));
               curNode = oss.str();
               oss.str("");

               // Add edge
               this->createAttr(pEdge, "color", color.at(i));
               oss << pathIt->edge(i).opId();
               this->createAttr(pEdge, "label", oss.str());
               oss.str("");
               pGraph->append_node(pEdge);

               if(pathIt->edge(i).arithId() == Arithmetics::SET)
               {
                  oss << "=";
               } else if(pathIt->edge(i).arithId() == Arithmetics::SETNEG)
               {
                  oss << "=-";
               } else if(pathIt->edge(i).arithId() == Arithmetics::ADD)
               {
                  oss << "+";
               } else if(pathIt->edge(i).arithId() == Arithmetics::SUB)
               {
                  oss << "-";
               }
               if(pathIt->edge(i).outId().at(0) != -1)
               {
                  if(pathIt->fieldId() == FieldType::SCALAR)
                  {
                     oss << std::endl << "Scalar";
                  } else if(pathIt->fieldId() == FieldType::CURL)
                  {
                     oss << std::endl << "Curl";
                  } else if(pathIt->fieldId() == FieldType::GRADIENT)
                  {
                     oss << std::endl << "Gradient";
                  } else if(pathIt->fieldId() == FieldType::GRADIENT2)
                  {
                     oss << std::endl << "Gradient2";
                  } else
                  {
                     oss << std::endl;
                  }

                  for(std::vector<int>::const_iterator outIt = pathIt->edge(i).outId().begin(); outIt != pathIt->edge(i).outId().end(); ++outIt)
                  {
                     if(dir == TransformDirection::FORWARD)
                     {
                        oss << " " << IoTools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(*outIt));
                     } else
                     {
                        oss << " " << IoTools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(*outIt));
                     }
                  }
               }
               this->createAttr(pNode, "label", oss.str());
               oss.str("");

               pGraph->append_node(pNode);
            }

            unid++;
         }
      }
   }

   void GxlWriter::graphTransformTree(const std::vector<Transform::TransformTree>& trees, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      // Set subgraph colors
      std::vector<std::string> color;
      color.push_back("blue");
      color.push_back("green");
      color.push_back("red");
      std::vector<std::string>::const_iterator colorIt = color.begin();

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the trees
      for(std::vector<Transform::TransformTree>::const_iterator treeIt = trees.begin(); treeIt != trees.end(); ++treeIt)
      {
         rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "t" << "_" << treeIt->name() << "_" << treeIt->comp<int>();
         std::string root = oss.str();
         oss.str("");
         pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(root.c_str(),0)));

         oss << IoTools::IdToHuman::toString(treeIt->name()) << " ";
         if(dir == TransformDirection::FORWARD)
         {
            oss << IoTools::IdToHuman::toString(treeIt->comp<FieldComponents::Physical::Id>());
         } else
         {
            oss << IoTools::IdToHuman::toString(treeIt->comp<FieldComponents::Spectral::Id>());
         }
         this->createAttr(pNode, "label", oss.str());
         oss.str("");

         pGraph->append_node(pNode);

         graphTransformTreeEdge(treeIt->root(), root, colorIt, pGraph, dir);
      }
   }

   void GxlWriter::graphTransformTreeEdge(const Transform::TransformTreeEdge& edge, const std::string& root, std::vector<std::string>::const_iterator colorIt, rapidxml::xml_node<> * pGraph, const TransformDirection::Id dir)
   {
      std::stringstream oss;

      for(Transform::TransformTreeEdge::EdgeType_citerator edgeIt = edge.edgeRange().first; edgeIt != edge.edgeRange().second; ++edgeIt)
      {
         rapidxml::xml_node<> *pNode = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << root << edgeIt->opId<int>();
         std::string nextRoot = oss.str();
         oss.str("");
         pNode->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(nextRoot.c_str(),0)));

         // Set "from" attribute
         rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
         pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(root.c_str(),0)));

         // Set "to" attribute
         pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(nextRoot.c_str(),0)));

         // Add edge
         this->createAttr(pEdge, "color", *colorIt);
         oss << edgeIt->opId<int>();
         if(edgeIt->recoverInput() && edgeIt->holdInput())
         {
            oss << std::endl << "(R,H)";
         } else if(edgeIt->recoverInput())
         {
            oss << std::endl << "(R)";
         } else if(edgeIt->holdInput())
         {
            oss << std::endl << "(H)";
         }
         if(edgeIt->recoverOutId() >= 0 && edgeIt->combinedOutId() >= 0)
         {
            oss << std::endl << "}" << edgeIt->recoverOutId() << "," << edgeIt->combinedOutId() << "{";
         } else if(edgeIt->recoverOutId() >= 0)
         {
            oss << std::endl << "}" << edgeIt->recoverOutId() << "{";
         } else if(edgeIt->combinedOutId() >= 0)
         {
            oss << std::endl << "}" << edgeIt->combinedOutId() << "{";
         }
         this->createAttr(pEdge, "label", oss.str());
         oss.str("");
         pGraph->append_node(pEdge);

         if(edgeIt->arithId() == Arithmetics::SET)
         {
            oss << "=";
         } else if(edgeIt->arithId() == Arithmetics::SETNEG)
         {
            oss << "=-";
         } else if(edgeIt->arithId() == Arithmetics::ADD)
         {
            oss << "+";
         } else if(edgeIt->arithId() == Arithmetics::SUB)
         {
            oss << "-";
         }
         if(edgeIt->outId<int>() != -1)
         {
            oss << std::endl;
            if(edgeIt->fieldId() == FieldType::SCALAR)
            {
               oss << "Scalar";
            } else if(edgeIt->fieldId() == FieldType::VECTOR)
            {
               // Just use name
            } else if(edgeIt->fieldId() == FieldType::CURL)
            {
               oss << "Curl";
            } else if(edgeIt->fieldId() == FieldType::GRADIENT)
            {
               oss << "Gradient";
            } else if(edgeIt->fieldId() == FieldType::GRADIENT2)
            {
               oss << "Gradient2";
            } else
            {
               throw Exception("Unknown field type ID in tree");
            }

            for(std::vector<int>::const_iterator outIt = edgeIt->outIds().begin(); outIt != edgeIt->outIds().end(); ++outIt)
            {
               if(dir == TransformDirection::FORWARD)
               {
                  oss << " " << IoTools::IdToHuman::toString(static_cast<FieldComponents::Spectral::Id>(*outIt));
               } else
               {
                  oss << " " << IoTools::IdToHuman::toString(static_cast<FieldComponents::Physical::Id>(*outIt));
               }
            }
         }
         this->createAttr(pNode, "label", oss.str());
         oss.str("");

         pGraph->append_node(pNode);

         colorIt++;
         graphTransformTreeEdge(*edgeIt, nextRoot, colorIt, pGraph, dir);
         colorIt--;
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
