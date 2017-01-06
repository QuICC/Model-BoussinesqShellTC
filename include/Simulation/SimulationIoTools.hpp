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
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "IoStats/IStatisticsAsciiEWriter.hpp"
#include "TypeSelectors/TransformCommSelector.hpp"

namespace QuICC {

   /**
    * @brief Implementation of the tools for IO related calculations
    */
   class SimulationIoTools
   {
      public:
         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<IoVariable::SharedIVariableAsciiEWriter>::iterator ascii_iterator;

         /// Typedef for an iterator over all the ASCII writers
         typedef std::vector<IoStats::SharedIStatisticsAsciiEWriter>::iterator stats_iterator;

         /**
          * @brief Update heavy calculations for ASCII file output
          */
         static void updateHeavyAscii(ascii_iterator asciiBegin, ascii_iterator asciiEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update pre calculations for statistics file output
          */
         static void updateStatsPre(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update calculations for statistics file output
          */
         static void updateStats(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

         /**
          * @brief Update post calculations for statistics file output
          */
         static void updateStatsPost(stats_iterator statsBegin, stats_iterator statsEnd, Transform::TransformCoordinatorType& coord);

      protected:

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
