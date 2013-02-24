/** \file ConfigurationWriter.cpp
 *  \brief Source of the implementation of the XML parameters file reader
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
#include "IoConfig/ConfigurationWriter.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string ConfigurationWriter::NAME = "parameters_template";

   ConfigurationWriter::ConfigurationWriter(const int dim, const std::string& type)
      : IConfigurationFile<XmlWriter>(dim, ConfigurationWriter::NAME, type)
   {
   }

   void ConfigurationWriter::write()
   {
      // Do pre writing processing
      this->preWrite();


      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         //
         // Write framework configuration
         //

         // Create framework XML node
         rapidxml::xml_node<> *pMaster = this->mXML.allocate_node(rapidxml::node_element, this->frameworkTag().c_str());
         this->mXML.append_node(pMaster);

         // Check if master node exists
         if(pMaster)
         {
            // Define component node pointer
            rapidxml::xml_node<> *pComponent;

            // Create master iterator
            std::map<FrameworkBlocks::Block, SharedConfigurationPart>::const_iterator itM;
            // Iterate over all master components
            for(itM = this->mFramework.begin(); itM != this->mFramework.end(); itM++)
            {
               // Create component node
               pComponent = this->mXML.allocate_node(rapidxml::node_element, itM->second->parent().c_str());
               pMaster->append_node(pComponent);

               // Check if component node exists
               if(pComponent)
               {
                  // Create integer component iterator
                  std::map<std::string, int>::const_iterator itIC;
                  // Iterate over all component entries
                  for(itIC = itM->second->iMap().begin(); itIC != itM->second->iMap().end(); itIC++)
                  {
                     // Create entry value
                     this->writeValue(itIC->second, pComponent, itIC->first);
                  }

                  // Create float component iterator
                  std::map<std::string, MHDFloat>::const_iterator itIF;
                  // Iterate over all component entries
                  for(itIF = itM->second->fMap().begin(); itIF != itM->second->fMap().end(); itIF++)
                  {
                     // Create entry value
                     this->writeValue(itIF->second, pComponent, itIF->first);
                  }

                  // Check if component data is correct
                  itM->second->checkData();
               } else
               {
                  throw Exception("Couldn't find " + itM->second->parent() + " component tag!");
               }
            }
         } else
         {
            throw Exception("Couldn't find " + this->frameworkTag() + " master tag!");
         }

         //
         // Write simulation configuration
         //
         if(this->mSimulation.size() > 0)
         {
            // Create framework XML node
            pMaster = this->mXML.allocate_node(rapidxml::node_element, this->simulationTag().c_str());
            this->mXML.append_node(pMaster);

            // Check if master node exists
            if(pMaster)
            {
               // Define component node pointer
               rapidxml::xml_node<> *pComponent;

               // Create master iterator
               std::map<SimulationBlocks::Block, SharedConfigurationPart>::const_iterator itM;
               // Iterate over all master components
               for(itM = this->mSimulation.begin(); itM != this->mSimulation.end(); itM++)
               {
                  // Create component node
                  pComponent = this->mXML.allocate_node(rapidxml::node_element, itM->second->parent().c_str());
                  pMaster->append_node(pComponent);

                  // Check if component node exists
                  if(pComponent)
                  {
                     // Create integer component iterator
                     std::map<std::string, int>::const_iterator itIC;
                     // Iterate over all component entries
                     for(itIC = itM->second->iMap().begin(); itIC != itM->second->iMap().end(); itIC++)
                     {
                        // Create entry value
                        this->writeValue(itIC->second, pComponent, itIC->first);
                     }

                     // Create float component iterator
                     std::map<std::string, MHDFloat>::const_iterator itIF;
                     // Iterate over all component entries
                     for(itIF = itM->second->fMap().begin(); itIF != itM->second->fMap().end(); itIF++)
                     {
                        // Create entry value
                        this->writeValue(itIF->second, pComponent, itIF->first);
                     }

                     // Check if component data is correct
                     itM->second->checkData();
                  } else
                  {
                     throw Exception("Couldn't find " + itM->second->parent() + " component tag!");
                  }
               }
            } else
            {
               throw Exception("Couldn't find " + this->simulationTag() + " master tag!");
            }
         }
      }

      // Write xml content to file
      this->mFile << this->mXML;

      // Share Parameters with other CPUs (if applicable)
      this->spreadParameters();

      // Do post writing processing
      this->postWrite();
   }
}
}
