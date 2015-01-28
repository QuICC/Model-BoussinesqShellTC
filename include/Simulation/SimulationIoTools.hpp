/** 
 * @file SimulationIoTools.hpp
 * @brief Implementation of the tools IO related calculations
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SIMULATIONIOTOOLS_HPP
#define SIMULATIONIOTOOLS_HPP

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/IVariableHeavyAsciiEWriter.hpp"
#include "TypeSelectors/TransformSelector.hpp"

namespace GeoMHDiSCC {

   /**
    * @brief Implementation of the tools for IO related calculations
    */
   class SimulationIoTools
   {
      public:
         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<IoVariable::SharedIVariableAsciiEWriter>::iterator ascii_iterator;

         /**
          * @brief Update heavy calculations for ASCII file output
          */
         static void updateHeavyAscii(ascii_iterator asciiBegin, ascii_iterator asciiEnd, Transform::TransformCoordinatorType& coord);

      protected:
         /**
          * @brief Compute heavy calculation for file
          */
         static void updateHeavyFile(SharedPtrMacro<IoVariable::IVariableHeavyAsciiEWriter> spAscii, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Do nothing in not heavy type
          */
         static void updateHeavyFile(SharedPtrMacro<IoVariable::IVariableAsciiEWriter> spAscii, Transform::TransformCoordinatorType& coord) {};


      private:
         /**
          * @brief Empty constructor
          */
         SimulationIoTools();

         /**
          * @brief Empty destructor
          */
         ~SimulationIoTools();
   };

}

#endif // SIMULATIONIOTOOLS_HPP
