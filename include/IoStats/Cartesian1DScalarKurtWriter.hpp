/** 
 * @file Cartesian1DScalarKurtWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) kurtosis calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSCALARKURTWRITER_HPP
#define CARTESIAN1DSCALARKURTWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoStats/IStatisticsAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"
#include "IoStats/Cartesian1DScalarAvgWriter.hpp"
#include "IoStats/Cartesian1DScalarRMSWriter.hpp"

namespace QuICC {

   namespace IoStats {

      /**
       * @brief Implementation of the ASCII Cartesian 1D (double periodic) skew calculation for a scalar field
       */
      class Cartesian1DScalarKurtWriter: public IStatisticsAsciiEWriter
      {
         public:
            /**
             * @brief Constructor
             *
             * @param prefix Prefix to use for file name
             * @param type Type of the file (typically scheme name)
             */
            Cartesian1DScalarKurtWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type);

            /**
             * @brief Destructor
             */
            virtual ~Cartesian1DScalarKurtWriter();

            /**
             * @brief Initialise the operator, transform and file
             */
            virtual void init();

            /**
             * @brief Compute kurtosis for scalar field
             */
            void compute(Transform::TransformCoordinatorType& coord);

            /**
             * @brief Post Compute for kurtosis of scalar field
             */
            void postCompute(Transform::TransformCoordinatorType& coord);

            /**
             * @brief Write State to file
             */
            virtual void write();

            /**
             * @brief Requires heavy calculation?
             */
            virtual bool isHeavy() const; 

         protected:

         private:
            /*
             *             * @brief Cartesian box area to normalize statistics
             *              */
            MHDFloat mArea;

            /**
             * @brief Storage for the scalar energy
             */
            SharedCartesian1DScalarAvgWriter mAvg;
            SharedCartesian1DScalarRMSWriter mRMS;
            Array mKurt;
      };

      inline bool Cartesian1DScalarKurtWriter::isHeavy() const
      {
         return true;
      }

      /// Typedef for a shared pointer of a HDF5 state file writer
      typedef SharedPtrMacro<Cartesian1DScalarKurtWriter> SharedCartesian1DScalarKurtWriter;

   }
}

#endif // CARTESIAN1DSCALARKURTWRITER_HPP
