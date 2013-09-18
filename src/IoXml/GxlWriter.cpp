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

   void GxlWriter::graphResolution(SharedResolution spRes)
   {
      if(spRes->cpu(0)->nDim() == 1)
      {
         this->graph1DResolution(spRes);

      } else if(spRes->cpu(0)->nDim() == 2)
      {
         this->graph1DResolution(spRes);

      } else if(spRes->cpu(0)->nDim() == 3)
      {
         this->graph3DResolution(spRes);
      }
   }

   void GxlWriter::graph1DResolution(SharedResolution spRes)
   {
      std::stringstream oss;

      // Set graph color
      std::string color = "blue";

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the CPUs
      for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
      {
         rapidxml::xml_node<> *pCpu = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "cpu" << cpu;
         pCpu->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
         oss.str("");

         oss << "CPU " << cpu;
         this->createAttr(pCpu, "label", oss.str());
         oss.str("");

         pGraph->append_node(pCpu);
      }

//      for(int j = 0; j < spRes->cpu(0)->nDim(); j++)
//      {
//      // Loop over CPUs
//      for(int i = 0; i < spRes->nCpu(); i++)
//      {
//         for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATF1D>(); m++)
//         {
//            stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATF1D>(m) << std::endl;
//         }
//      }
//      stream << "- -" << std::endl;
//
//      for(int i = 0; i < spRes->nCpu(); i++)
//      {
//         for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATB1D>(); m++)
//         {
//            stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATB1D>(m) << std::endl;
//         }
//      }
//      }
   }

   void GxlWriter::graph2DResolution(SharedResolution spRes)
   {
      std::stringstream oss;

      // Set graph color
      std::string color = "blue";

      // Get master GXL tag
      rapidxml::xml_node<> *pGxl = this->mXML.first_node(this->GXLTAG.c_str());

      // Create master graph
      rapidxml::xml_node<> *pGraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
      pGraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string("master",0)));
      pGxl->append_node(pGraph);

      // Loop over the CPUs
      for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
      {
         rapidxml::xml_node<> *pCpu = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "cpu" << cpu;
         pCpu->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
         oss.str("");

         oss << "CPU " << cpu;
         this->createAttr(pCpu, "label", oss.str());
         oss.str("");

         pGraph->append_node(pCpu);
      }
//      for(int j = 0; j < spRes->cpu(0)->nDim(); j++)
//      {
//         // Loop over CPUs
//         for(int i = 0; i < spRes->nCpu(); i++)
//         {
//            for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(); l++)
//            {
//               for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATF1D>(l); m++)
//               {
//                  stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATF1D>(m,l) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l) << std::endl;
//               }
//            }
//         }
//         stream << "- - -" << std::endl;
//
//         for(int i = 0; i < spRes->nCpu(); i++)
//         {
//            for(int l = 0; l < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DAT2D>(); l++)
//            {
//               for(int m = 0; m < spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->dim<Dimensions::Data::DATB1D>(l); m++)
//               {
//                  stream << i << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DATB1D>(m,l) << " " << spRes->cpu(i)->dim(static_cast<Dimensions::Transform::Id>(j))->idx<Dimensions::Data::DAT2D>(l) << std::endl;
//               }
//            }
//         }
//         stream << "- - -" << std::endl;
//      }
   }

   void GxlWriter::graph3DResolution(SharedResolution spRes)
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
      for(int i = 0; i < 2; i++)
      {
         // Cast integer to dimension ID
         Dimensions::Transform::Id fwdDim = static_cast<Dimensions::Transform::Id>(i);

         // Create forward coordinates list
         std::vector<std::set<std::tr1::tuple<int,int,int> > > coordFwd;
         coordFwd.reserve(spRes->nCpu());

         // Loop over all CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Initialise coordinate list for current CPU
            coordFwd.push_back(std::set<std::tr1::tuple<int,int,int> >());

            // Loop over slow data dimension
            int i_, j_, k_;
            for(int k = 0; k < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Get global slow index
               k_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DAT3D>(k);

               // Loop over middle data dimension
               for(int j = 0; j < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DAT2D>(k); ++j)
               {
                  // Get global middle index
                  j_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DAT2D>(j,k);

                  // Loop over middle Forward data dimension
                  for(int i = 0; i < spRes->cpu(cpu)->dim(fwdDim)->dim<Dimensions::Data::DATF1D>(k); ++i)
                  {
                     // Get global fast index
                     i_ = spRes->cpu(cpu)->dim(fwdDim)->idx<Dimensions::Data::DATF1D>(i,k);

                     // Add coordinate to list
                     coordFwd.back().insert(std::tr1::make_tuple(i_,k_,j_));
                  }
               }
            }
         }

         // Create backward coordinate list
         std::vector<std::set<std::tr1::tuple<int,int,int> > > coordBwd;
         coordBwd.reserve(spRes->nCpu());

         // Loop over all CPUs
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            // Initialise coordinate list for current CPU
            coordBwd.push_back(std::set<std::tr1::tuple<int,int,int> >());

            // Loop over slow data dimension
            int i_, j_, k_;
            for(int k = 0; k < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT3D>(); ++k)
            {
               // Get global slow index
               k_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT3D>(k);

               // Loop over middle backward data dimension
               for(int i = 0; i < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DATB1D>(k); ++i)
               {
                  // Get global fast index
                  i_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DATB1D>(i,k);

                  // Loop over middle data dimension
                  for(int j = 0; j < spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->dim<Dimensions::Data::DAT2D>(k); ++j)
                  {
                     // Get global middle index
                     j_ = spRes->cpu(cpu)->dim(Dimensions::jump(fwdDim,1))->idx<Dimensions::Data::DAT2D>(j,k);

                     // Add coordinate to list
                     coordBwd.back().insert(coordBwd.back().end()--,std::tr1::make_tuple(j_,i_,k_));
                  }
               }
            }
         }

         // Create subgraph
         rapidxml::xml_node<> *pNSubgraph = this->mXML.allocate_node(rapidxml::node_element, "node");
         oss << "N_" << lName.at(i);
         pNSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
         oss.str("");
         rapidxml::xml_node<> *pSubgraph = this->mXML.allocate_node(rapidxml::node_element, "graph");
         pSubgraph->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(lName.at(i).c_str(),0)));

         if(fwdDim == Dimensions::Transform::TRA1D)
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
         for(int cpu = 0; cpu < spRes->nCpu(); cpu++)
         {
            rapidxml::xml_node<> *pCpu = this->mXML.allocate_node(rapidxml::node_element, "node");
            oss << "cpu" << cpu << sName.at(i);
            pCpu->append_attribute(this->mXML.allocate_attribute("id", this->mXML.allocate_string(oss.str().c_str(),0)));
            oss.str("");

            oss << "CPU " << cpu;
            this->createAttr(pCpu, "label", oss.str());
            oss.str("");

            pSubgraph->append_node(pCpu);
         }

         // Create set for shared coordinates
         std::set<std::tr1::tuple<int,int,int> > shared;

         // Loop over forward cpus
         for(int cpu1D = 0; cpu1D < spRes->nCpu(); cpu1D++)
         {
            // Loop over backward cpus
            for(int cpu2D = 0; cpu2D < spRes->nCpu(); cpu2D++)
            {
               // Compute coordinate intersection
               std::set_intersection(coordFwd.at(cpu1D).begin(), coordFwd.at(cpu1D).end(), coordBwd.at(cpu2D).begin(), coordBwd.at(cpu2D).end(), std::inserter(shared, shared.begin()));

               // Create link if at least one coordinate is shared
               if(shared.size() > 0)
               {
                  rapidxml::xml_node<> *pEdge = this->mXML.allocate_node(rapidxml::node_element, "edge");
                  oss << "cpu" << cpu1D << sName.at(i);
                  pEdge->append_attribute(this->mXML.allocate_attribute("from", this->mXML.allocate_string(oss.str().c_str(),0)));
                  oss.str("");
                  oss << "cpu" << cpu2D << sName.at(i);
                  pEdge->append_attribute(this->mXML.allocate_attribute("to", this->mXML.allocate_string(oss.str().c_str(),0)));
                  oss.str("");

                  this->createAttr(pEdge, "color", color.at(i));

                  pSubgraph->append_node(pEdge);
               }

               // Clear shared list
               shared.clear();
            }
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
