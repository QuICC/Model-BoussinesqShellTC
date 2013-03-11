/** \file ConfigurationReader.cpp
 *  \brief Source of the implementation of the XML parameters file reader
 */

// Configuration includes
//
#include "Framework/FrameworkMacro.h"

// System includes
//

// External includes
//

// Class include
//
#include "IoConfig/ConfigurationReader.hpp"

// Project includes
//
#include "Exceptions/Exception.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   const std::string ConfigurationReader::NAME = "parameters";

   ConfigurationReader::ConfigurationReader(const int dim, const std::string& type)
      : IConfigurationFile<IoXml::IXmlReader>(dim, ConfigurationReader::NAME, type)
   {
   }

   ConfigurationReader::~ConfigurationReader()
   {
   }

   void ConfigurationReader::read()
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         //
         // Read framework configuration
         //

         // Get master pointer to framework XML code
         rapidxml::xml_node<> *pMaster = this->mXML.first_node(this->frameworkTag().c_str());

         std::pair<std::map<std::string,int>::const_iterator,std::map<std::string,int>::const_iterator> iRange;
         std::pair<std::map<std::string,MHDFloat>::const_iterator,std::map<std::string,MHDFloat>::const_iterator> fRange;

         // Check if master node exists
         if(pMaster)
         {
            // Define component node pointer
            rapidxml::xml_node<> *pComponent;

            // Create master iterator
            std::map<FrameworkBlocks::Id, SharedIConfigurationPart>::iterator itM;
            // Iterate over all master components
            for(itM = this->mFramework.begin(); itM != this->mFramework.end(); itM++)
            {
               // Get pointer to component node
               pComponent = pMaster->first_node(itM->second->parent().c_str());

               // Check if component node exists
               if(pComponent)
               {
                  // Create integer component iterator
                  std::map<std::string, int>::const_iterator itIC;
                  iRange = itM->second->iRange();
                  // Iterate over all component entries
                  for(itIC = iRange.first; itIC != iRange.second; itIC++)
                  {
                     // Read entry value from file
                     int val;
                     this->readValue(val, pComponent, itIC->first);

                     // Store value
                     itM->second->setIValue(itIC->first, val);
                  }

                  // Create float component iterator
                  std::map<std::string, MHDFloat>::const_iterator itIF;
                  fRange = itM->second->fRange();
                  // Iterate over all component entries
                  for(itIF = fRange.first; itIF != fRange.second; itIF++)
                  {
                     // Read entry value from file
                     MHDFloat val;
                     this->readValue(val, pComponent, itIF->first);

                     // Store value
                     itM->second->setFValue(itIF->first, val);
                  }

                  // Check if component data is correct
                  itM->second->checkData();
               } else
               {
                  throw Exception("Couldn't find " + itM->second->parent() + " component tag from configuration file!");
               }
            }
         } else
         {
            throw Exception("Couldn't find " + this->frameworkTag() + " master tag from configuration file!");
         }

         //
         // Read simulation configuration
         //
         if(this->mSimulation.size() > 0)
         {
            // Get master pointer to framework XML code
            pMaster = this->mXML.first_node(this->simulationTag().c_str());

            // Check if master node exists
            if(pMaster)
            {
               // Define component node pointer
               rapidxml::xml_node<> *pComponent;

               // Create master iterator
               std::map<SimulationBlocks::Id, SharedIConfigurationPart>::iterator itM;
               // Iterate over all master components
               for(itM = this->mSimulation.begin(); itM != this->mSimulation.end(); itM++)
               {
                  // Get pointer to component node
                  pComponent = pMaster->first_node(itM->second->parent().c_str());

                  // Check if component node exists
                  if(pComponent)
                  {
                     // Create integer component iterator
                     std::map<std::string, int>::const_iterator itIC;
                     iRange = itM->second->iRange();
                     // Iterate over all component entries
                     for(itIC = iRange.first; itIC != iRange.second; itIC++)
                     {
                        // Read entry value from file
                        int value;
                        this->readValue(value, pComponent, itIC->first);

                        // Store value
                        itM->second->setIValue(itIC->first, value);
                     }

                     // Create float component iterator
                     std::map<std::string, MHDFloat>::const_iterator itIF;
                     fRange = itM->second->fRange();
                     // Iterate over all component entries
                     for(itIF = fRange.first; itIF != fRange.second; itIF++)
                     {
                        // Read entry value from file
                        MHDFloat value;
                        this->readValue(value, pComponent, itIF->first);

                        // Store value
                        itM->second->setFValue(itIF->first, value);
                     }

                     // Check if component data is correct
                     itM->second->checkData();
                  } else
                  {
                     throw Exception("Couldn't find " + itM->second->parent() + " component tag from configuration file!");
                  }
               }
            } else
            {
               throw Exception("Couldn't find " + this->simulationTag() + " master tag from configuration file!");
            }
         }
      }

      // Share Parameters with other CPUs (if applicable)
      this->spreadParameters();
   }
}
}
