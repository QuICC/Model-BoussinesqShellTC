/** \file IConfigurationFile.hpp 
 *  \brief Implementation of the base for a configuration file
 *
 *  \mhdBug Needs test
 */

#ifndef ICONFIGURATIONFILE_HPP
#define ICONFIGURATIONFILE_HPP

// Configuration includes
//
#include "Framework/FrameworkMacro.h"
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//
#include <string>
#include <map>

// External includes
//

// Project includes
//
#include "Exceptions/Exception.hpp"
#include "Base/Typedefs.hpp"
#include "IoTools/Formatter.hpp"
#include "IoConfig/ConfigParts/IConfigurationPart.hpp"
#include "IoConfig/ConfigParts/TruncationPart.hpp"
#include "IoConfig/ConfigParts/ParallelPart.hpp"
#include "IoConfig/ConfigParts/TimesteppingPart.hpp"
#include "IoConfig/ConfigParts/RunPart.hpp"
#include "IoConfig/ConfigParts/IoPart.hpp"

namespace GeoMHDiSCC {

namespace IoConfig {

   /**
    * @brief Simple struck holding the framework configuration blocks
    */
   struct FrameworkBlocks {

      /**
       * @name Enum for framework configuration block
       */
      enum Id {TRUNCATION, PARALLEL, TIMESTEPPING, RUN, IO};
   };

   /**
    * @brief Simple struck holding the simulation configuration blocks
    */
   struct SimulationBlocks {

      /**
       * @name Enum of the possible blocks
       */
      enum Id {PHYSICAL, BOUNDARY};
   };

   /**
    * @brief Implementation of the base for a configuration file
    */
   template <typename TBase> class IConfigurationFile: public TBase
   {
      public:
         /**
          * @brief Constructor
          *
          * @param dim  Dimensionality of simulation
          * @param name Name of the file
          * @param type Type of the simulation
          */
         IConfigurationFile(const int dim, const std::string& name, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~IConfigurationFile();

         /**
          * @brief Output run information
          */
         void printInfo() const;

         /**
          * @brief Add configuration part to simulation configuration
          *
          * @param id      ID of the simulation configuration block
          * @param spPart  Configuration part to add
          */
         void addPart(SimulationBlocks::Id id, SharedIConfigurationPart spPart);

         /**
          * @brief Get truncation part
          */
         SharedCIConfigurationPart spTruncation() const;

         /**
          * @brief Get parallel part
          */
         SharedCIConfigurationPart spParallel() const;

         /**
          * @brief Get timestepping part
          */
         SharedCIConfigurationPart spTimestepping() const;

         /**
          * @brief Get run part
          */
         SharedCIConfigurationPart spRun() const;

         /**
          * @brief Get run part
          */
         SharedCIConfigurationPart spIo() const;

         /**
          * @brief Get physical part
          */
         SharedCIConfigurationPart spPhysical() const;

         /**
          * @brief Get boundary part
          */
         SharedCIConfigurationPart spBoundary() const;
         
      protected:
         /**
          * @brief Get name of the framework XML tag
          */
         std::string frameworkTag() const;

         /**
          * @brief Get name of the framework XML tag
          */
         std::string simulationTag() const;

         /**
          * @brief Framework part of configuration file
          */
         std::map<FrameworkBlocks::Id, SharedIConfigurationPart>  mFramework;

         /**
          * @brief Simulation part of configuration file
          */
         std::map<SimulationBlocks::Id, SharedIConfigurationPart>  mSimulation;

         /**
          * @brief Initialise the different parts of the framework part
          *
          * @param dim Dimensionality of simulation
          */
         void initFramework(const int dim);

         /**
          * @brief Spread parameters over parallel simulation
          *
          * \mhdBug MPI code does not handle different MHDFloat types, only double
          */
         void spreadParameters();

      private:
         /**
          * @brief Name of the framework XML tag
          */
         static const std::string FRAMEWORKTAG;

         /**
          * @brief Name of the simulation XML tag
          */
         static const std::string SIMULATIONTAG;

         /**
          * @brief EXTENSION of the configuration file
          */
         static const std::string EXTENSION;

         /**
          * @brief HEADER of the configuration file
          */
         static const std::string HEADER;

         /**
          * @brief VERSION of the configuration file
          */
         static const std::string VERSION;
   };

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spTruncation() const
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::TRUNCATION) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::TRUNCATION)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spParallel() const
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::PARALLEL) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::PARALLEL)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spTimestepping() const
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::TIMESTEPPING) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::TIMESTEPPING)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spRun() const
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::RUN) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::RUN)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spIo() const
   {
      // Make sure initialisation was correct
      assert(this->mFramework.find(FrameworkBlocks::IO) != this->mFramework.end());

      return this->mFramework.find(FrameworkBlocks::IO)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spPhysical() const
   {
      // Make sure initialisation was correct
      assert(this->mSimulation.find(SimulationBlocks::PHYSICAL) != this->mSimulation.end());

      return this->mSimulation.find(SimulationBlocks::PHYSICAL)->second;
   }

   template <typename TBase> inline SharedCIConfigurationPart IConfigurationFile<TBase>::spBoundary() const
   {
      // Make sure initialisation was correct
      assert(this->mSimulation.find(SimulationBlocks::BOUNDARY) != this->mSimulation.end());

      return this->mSimulation.find(SimulationBlocks::BOUNDARY)->second;
   }

   template <typename TBase> std::string IConfigurationFile<TBase>::frameworkTag() const
   {
      return IConfigurationFile<TBase>::FRAMEWORKTAG;
   }

   template <typename TBase> std::string IConfigurationFile<TBase>::simulationTag() const
   {
      return IConfigurationFile<TBase>::SIMULATIONTAG;
   }

   template <typename TBase> const std::string IConfigurationFile<TBase>::FRAMEWORKTAG = "framework";

   template <typename TBase> const std::string IConfigurationFile<TBase>::SIMULATIONTAG = "simulation";

   template <typename TBase> const std::string IConfigurationFile<TBase>::EXTENSION = ".cfg";

   template <typename TBase> const std::string IConfigurationFile<TBase>::HEADER = "ConfigurationFile";

   template <typename TBase> const std::string IConfigurationFile<TBase>::VERSION = "1.0";

   template <typename TBase> IConfigurationFile<TBase>::IConfigurationFile(const int dim, const std::string& name, const std::string& type)
      : TBase(name, IConfigurationFile<TBase>::EXTENSION, IConfigurationFile<TBase>::HEADER, type, IConfigurationFile<TBase>::VERSION)
   {
      // Initialise the file descriptors
      this->initFramework(dim); 
   }

   template <typename TBase> IConfigurationFile<TBase>::~IConfigurationFile()
   {
   }

   template <typename TBase> void IConfigurationFile<TBase>::initFramework(const int dim)
   {
      // Create shared pointer
      SharedIConfigurationPart spPart;

      //
      // Create framework content
      
      // Add truncation part
      spPart = SharedTruncationPart(new TruncationPart(dim));
      this->mFramework.insert(std::make_pair(FrameworkBlocks::TRUNCATION, spPart));
      
      // Add parallel part
      spPart = SharedParallelPart(new ParallelPart());
      this->mFramework.insert(std::make_pair(FrameworkBlocks::PARALLEL, spPart));
      
      // Add parallel part
      spPart = SharedTimesteppingPart(new TimesteppingPart());
      this->mFramework.insert(std::make_pair(FrameworkBlocks::TIMESTEPPING, spPart));
      
      // Add parallel part
      spPart = SharedRunPart(new RunPart());
      this->mFramework.insert(std::make_pair(FrameworkBlocks::RUN, spPart));
      
      // Add IO part
      spPart = SharedIoPart(new IoPart());
      this->mFramework.insert(std::make_pair(FrameworkBlocks::IO, spPart));
   }

   template <typename TBase> void IConfigurationFile<TBase>::addPart(SimulationBlocks::Id id, SharedIConfigurationPart spPart)
   {
      this->mSimulation.insert(std::make_pair(id, spPart));
   }

   template <typename TBase> void IConfigurationFile<TBase>::printInfo() const
   {
      // Check if the framework allows IO to be performed
      if(FrameworkMacro::allowsIO())
      {
         // Create output header
         IoTools::Formatter::printNewline(std::cout);
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printCentered(std::cout, "Configuration parameters", '*');
         IoTools::Formatter::printLine(std::cout, '-');
         IoTools::Formatter::printNewline(std::cout);

         // Create framework constant part iterator
         std::map<FrameworkBlocks::Id, SharedIConfigurationPart>::const_iterator itF;
         // Iterate over all master components of framework configuration
         for(itF = this->mFramework.begin(); itF != this->mFramework.end(); itF++)
         {
            itF->second->printInfo();
         }

         // Create simulation constant part iterator
         std::map<SimulationBlocks::Id, SharedIConfigurationPart>::const_iterator itS;
         // Iterate over all master components of simulation configuration
         for(itS = this->mSimulation.begin(); itS != this->mSimulation.end(); itS++)
         {
            itS->second->printInfo();
         }
      }
   }

   template <typename TBase> void IConfigurationFile<TBase>::spreadParameters()
   {
      //
      // Start of MPI block
      //
      #ifdef GEOMHDISCC_MPI

      // Create MPI compatible storage for the integer values
      std::vector<int> iData;
      // Create MPI compatible storage for the float values
      std::vector<MHDFloat> fData;

      // 
      // Gather data
      //
      
      // Create framework iterator
      std::map<FrameworkBlocks::Id, SharedIConfigurationPart>::const_iterator citF;
      // Iterate over all master components
      for(citF = this->mFramework.begin(); citF != this->mFramework.end(); citF++)
      {
         // Create integer component iterator
         std::map<std::string, int>::const_iterator citIC;
         // Iterate over all component entries
         for(citIC = citF->second->iMap().begin(); citIC != citF->second->iMap().end(); citIC++)
         {
            iData.push_back(citIC->second);
         }

         // Create float component iterator
         std::map<std::string, MHDFloat>::const_iterator citIF;
         // Iterate over all component entries
         for(citIF = citF->second->fMap().begin(); citIF != citF->second->fMap().end(); citIF++)
         {
            fData.push_back(citIF->second);
         }
      }

      // Create simulation iterator
      std::map<SimulationBlocks::Id, SharedIConfigurationPart>::const_iterator citS;
      // Iterate over all master components of the simulation data
      for(citS = this->mSimulation.begin(); citS != this->mSimulation.end(); citS++)
      {
         // Create integer component iterator
         std::map<std::string, int>::const_iterator citIC;
         // Iterate over all component entries
         for(citIC = citS->second->iMap().begin(); citIC != citS->second->iMap().end(); citIC++)
         {
            iData.push_back(citIC->second);
         }

         // Create float component iterator
         std::map<std::string, MHDFloat>::const_iterator citIF;
         // Iterate over all component entries
         for(citIF = citS->second->fMap().begin(); citIF != citS->second->fMap().end(); citIF++)
         {
            fData.push_back(citIF->second);
         }
      }

      // Various MPI broadcast data
      int nBlocks = 2;
      MPI_Aint    displ[nBlocks];
      int         blocks[nBlocks];
      MPI_Datatype   types[nBlocks];
      MPI_Aint    element;

      // Create integer data part
      MPI_Get_address(&iData[0], &element);
      displ[0] = element;
      blocks[0] = iData.size();
      types[0] = MPI_INT;

      // Create float data part
      MPI_Get_address(&fData[0], &element);
      displ[1] = element;
      blocks[1] = fData.size();
      types[1] = MPI_DOUBLE;

      // MPI data type for the combined integer and float data
      MPI_Datatype   cfgType;

      // Create MPI datatype
      MPI_Type_create_struct(nBlocks, blocks, displ, types, &cfgType);
      // Commit MPI datatype
      MPI_Type_commit(&cfgType);

      // Broadcast the information
      MPI_Bcast(MPI_BOTTOM, 1, cfgType, FrameworkMacro::IO_RANK, MPI_COMM_WORLD);
      
      // Free the datatype
      MPI_Type_free(&cfgType);

      // 
      // Distribute data
      //
      
      // Global integer index
      int iIdx = 0;
      // Global float index
      int fIdx = 0;

      std::pair<std::map<std::string, int>::const_iterator,std::map<std::string, int>::const_iterator> iRange;
      std::pair<std::map<std::string, MHDFloat>::const_iterator,std::map<std::string, MHDFloat>::const_iterator> fRange;

      // Create framework iterator
      std::map<FrameworkBlocks::Id, SharedIConfigurationPart>::iterator itF;
      // Iterate over all master components
      for(itF = this->mFramework.begin(); itF != this->mFramework.end(); itF++)
      {
         // Create integer component iterator
         std::map<std::string, int>::const_iterator itIC;
         iRange = itF->second->iRange();
         // Iterate over all component entries
         for(itIC = iRange.first; itIC != iRange.second; itIC++)
         {
            itF->second->setIValue(itIC->first, iData.at(iIdx));
            // Increment integer index
            iIdx++;
         }

         // Create float component iterator
         std::map<std::string, MHDFloat>::const_iterator itIF;
         fRange = itF->second->fRange();
         // Iterate over all component entries
         for(itIF = fRange.first; itIF != fRange.second; itIF++)
         {
            itF->second->setFValue(itIF->first, fData.at(fIdx));
            // Increment float index
            fIdx++;
         }
      }

      // Create simulation iterator
      std::map<SimulationBlocks::Id, SharedIConfigurationPart>::iterator itS;
      // Iterate over all master components of the simulation data
      for(itS = this->mSimulation.begin(); itS != this->mSimulation.end(); itS++)
      {
         // Create integer component iterator
         std::map<std::string, int>::const_iterator itIC;
         iRange = itS->second->iRange();
         // Iterate over all component entries
         for(itIC = iRange.first; itIC != iRange.second; itIC++)
         {
            itS->second->setIValue(itIC->first, iData.at(iIdx));
            // Increment integer index
            iIdx++;
         }

         // Create float component iterator
         std::map<std::string, MHDFloat>::const_iterator itIF;
         fRange = itS->second->fRange();
         // Iterate over all component entries
         for(itIF = fRange.first; itIF != fRange.second; itIF++)
         {
            itS->second->setFValue(itIF->first, fData.at(fIdx));
            // Increment float index
            fIdx++;
         }
      }

      //
      // End of MPI block
      //
      #endif // GEOMHDISCC_MPI
   }

}
}

#endif // ICONFIGURATIONFILE_HPP
