/** 
 * @file Cartesian1DScalarSkewWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSCALARSKEWWRITER_HPP
#define CARTESIAN1DSCALARSKEWWRITER_HPP

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

namespace GeoMHDiSCC {

namespace IoStats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) skew calculation for a scalar field
    */
   class Cartesian1DScalarSkewWriter: public IStatisticsAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarSkewWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarSkewWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute skew for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Post Compute skew for scalar field
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

          MHDFloat mArea;
         /**
          * @brief Storage for the scalar energy
          */
          SharedCartesian1DScalarAvgWriter mAvg;
          SharedCartesian1DScalarRMSWriter mRMS;
         Array mSkew;
   };

   inline bool Cartesian1DScalarSkewWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DScalarSkewWriter> SharedCartesian1DScalarSkewWriter;

}
}

#endif // CARTESIAN1DSCALARSKEWWRITER_HPP
