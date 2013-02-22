/** \file ExampleModel.hpp
 *  \brief Implementation of an example physical model
 */

#ifndef EXAMPLEMODEL_HPP
#define EXAMPLEMODEL_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

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
          * @brief Create a shared simulation for the model
          */
         static SharedSimulation createSimulation();

      protected:
         /**
          * @brief Get the required equations
          */
         static void equations();

         /**
          * @brief Get the required XML input files
          */
         static void xmlInputs();

         /**
          * @brief Get the required HDF5 input files
          */
         static void hdf5Inputs();

         /**
          * @brief Get the required ASCII output files
          */
         static void asciiOutputs();

         /**
          * @brief Get the required HDF5 output files
          */
         static void hdf5Outputs();


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
