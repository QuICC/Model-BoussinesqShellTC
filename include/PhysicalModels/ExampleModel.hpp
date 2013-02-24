/** \file ExampleModel.hpp
 *  \brief Implementation of an example physical model
 */

#ifndef EXAMPLEMODEL_HPP
#define EXAMPLEMODEL_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Simulation/Simulation.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of an example physical model
    */
   class ExampleModel
   {
      public:
         /**
          * @brief Get the required equations
          *
          * @param spSim   Shared simulation object
          */
         static void addEquations(SharedSimulation spSim);

         /**
          * @brief Get the required XML input files
          *
          * @param spSim   Shared simulation object
          */
         static void setConfigurationFile(SharedSimulation spSim);

         /**
          * @brief Get the required HDF5 input files
          *
          * @param spSim   Shared simulation object
          */
         static void setInitialStateFile(SharedSimulation spSim);

         /**
          * @brief Get the required ASCII output files
          *
          * @param spSim   Shared simulation object
          */
         static void addAsciiOutputFiles(SharedSimulation spSim);

         /**
          * @brief Get the required HDF5 output files
          *
          * @param spSim   Shared simulation object
          */
         static void addHdf5OutputFiles(SharedSimulation spSim);

      protected:

      private:
         /**
          * @brief Constructor
          */
         ExampleModel();

         /**
          * @brief Simple empty destructor
          */
         ~ExampleModel();
   };

}

#endif // EXAMPLEMODEL_HPP
