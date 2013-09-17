/** 
 * @file ConfigurationWriter.cpp
 * @brief Source of the implementation of the XML parameters file reader
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
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

   const std::string ConfigurationWriter::NAME = "parameters_template_";

   ConfigurationWriter::ConfigurationWriter(const int dim, const std::vector<bool>& isPeriodicBox, const std::string& type)
      : IConfigurationFile<IoXml::IXmlWriter>(dim, isPeriodicBox, ConfigurationWriter::NAME+type, type)
   {
   }

   ConfigurationWriter::~ConfigurationWriter()
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

         std::pair<std::map<std::string,int>::const_iterator,std::map<std::string,int>::const_iterator> iRange;
         std::pair<std::map<std::string,MHDFloat>::const_iterator,std::map<std::string,MHDFloat>::const_iterator> fRange;

         // Check if master node exists
         if(pMaster)
         {
            // Define component node pointer
            rapidxml::xml_node<> *pComponent;

            // Create master iterator
            std::map<FrameworkBlocks::Id, SharedIConfigurationPart>::const_iterator itM;
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
                  iRange = itM->second->iRange();
                  // Iterate over all component entries
                  for(itIC = iRange.first; itIC != iRange.second; itIC++)
                  {
                     // Create entry value
                     this->writeValue(itIC->second, pComponent, itIC->first);
                  }

                  // Create float component iterator
                  std::map<std::string, MHDFloat>::const_iterator itIF;
                  fRange = itM->second->fRange();
                  // Iterate over all component entries
                  for(itIF = fRange.first; itIF != fRange.second; itIF++)
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
               std::map<SimulationBlocks::Id, SharedIConfigurationPart>::const_iterator itM;
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
                     iRange = itM->second->iRange();
                     // Iterate over all component entries
                     for(itIC = iRange.first; itIC != iRange.second; itIC++)
                     {
                        // Create entry value
                        this->writeValue(itIC->second, pComponent, itIC->first);
                     }

                     // Create float component iterator
                     std::map<std::string, MHDFloat>::const_iterator itIF;
                     fRange = itM->second->fRange();
                     // Iterate over all component entries
                     for(itIF = fRange.first; itIF != fRange.second; itIF++)
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

   SharedIConfigurationPart ConfigurationWriter::rspTruncation()
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::TRUNCATION) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::TRUNCATION)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspParallel()
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::PARALLEL) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::PARALLEL)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspTimestepping()
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::TIMESTEPPING) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::TIMESTEPPING)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspRun()
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::RUN) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::RUN)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspIo()
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::IO) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::IO)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspPhysical()
   {
      // Make sure initialisation was correct
      assert(this->mSimulation.find(SimulationBlocks::PHYSICAL) != this->mSimulation.end());

      return this->mSimulation.find(SimulationBlocks::PHYSICAL)->second;
   }

   SharedIConfigurationPart ConfigurationWriter::rspBoundary()
   {
      // Make sure initialisation was correct
      assert(this->mSimulation.find(SimulationBlocks::BOUNDARY) != this->mSimulation.end());

      return this->mSimulation.find(SimulationBlocks::BOUNDARY)->second;
   }
}
}
