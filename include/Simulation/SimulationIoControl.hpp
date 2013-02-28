/** \file SimulationIoControl.hpp
 *  \brief Implementation of the control instance for the simulation IO related operations
 *
 *  \mhdBug Needs test
 */

#ifndef SIMULATIONIOCONTROL_HPP
#define SIMULATIONIOCONTROL_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IoAscii/IAsciiWriter.hpp"
#include "IoAscii/StdOutPipe.hpp"
#include "IoHdf5/IHdf5Writer.hpp"
#include "IoConfig/ConfigurationReader.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the control instance for the simulation IO related operations
    */
   class SimulationIoControl
   {
      public:
         /**
          * @brief Constructor
          */
         SimulationIoControl();

         /**
          * @brief Destructor
          */
         ~SimulationIoControl();

         /**
          * @brief Set the configuration file
          */
         void setConfigurationFile(IoConfig::SharedConfigurationReader spCfgFile);

         /**
          * @brief Initialise the IO system
          */
         void init();

         /**
          * @brief Cleanup unused memory from IO system
          */
         void cleanup();

         /**
          * @brief Finalise the IO system
          */
         void finalize();

         /**
          * @brief Initialise the writer files
          */
         void initWriters();

         /**
          * @brief Finalise the writer files
          */
         void finalizeWriters();

         /**
          * @brief Add an ASCII output file
          *
          * @param spOutFile Shared ASCII writer
          */
         void addOutputFile(IoAscii::SharedIAsciiWriter spOutFile);

         /**
          * @brief Add a HDF5 output file
          *
          * @param spOutFile Shared HDF5 writer
          */
         void addOutputFile(IoHdf5::SharedIHdf5Writer  spOutFile);

         /**
          * @brief Write ASCII data
          */
         void writeAscii();

         /**
          * @brief Write HDF5 data
          */
         void writeHdf5();
         
      protected:

      private:
         /**
          * @brief Vector of ASCII output files
          */
         std::vector<IoAscii::SharedIAsciiWriter> mAsciiWriters;

         /**
          * @brief Vector of HDF5 output files
          */
         std::vector<IoHdf5::SharedIHdf5Writer> mHdf5Writers;

         /**
          * @brief Handle to StdMessage buffer
          */
         IoAscii::SharedStdOutPipe   mspStdOut;

         /**
          * @brief Handle to the configuration file
          */
         IoConfig::SharedConfigurationReader    mspCfgFile;

         /**
          * @brief Initialise the configuration file
          */
         void initCfg();

         /**
          * @brief Initialise the StdMessage file
          */
         void initStdOut();
   };

}

#endif // SIMULATIONIOCONTROL_HPP
