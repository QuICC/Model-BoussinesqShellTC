/** 
 * @file Cartesian1DScalarRMSWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSCALARRMSWRITER_HPP
#define CARTESIAN1DSCALARRMSWRITER_HPP

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

namespace QuICC {

namespace IoStats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
    */
   class Cartesian1DScalarRMSWriter: public IStatisticsAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarRMSWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarRMSWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute RMS for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief PostCompute RMS for scalar field
          */
         void postCompute(Transform::TransformCoordinatorType& coord);
         
         /**
          * @brief Write State to file
          */
         virtual void write();

         /**
          * To shrare RMS with other stats
          */
         const Array & RMS() const;

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:

      private:
         /*
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mArea;

         /**
          * @brief Storage for the scalar energy
          */
         SharedCartesian1DScalarAvgWriter mAvg;
         Array mRMS;
   };

   inline bool Cartesian1DScalarRMSWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DScalarRMSWriter> SharedCartesian1DScalarRMSWriter;

}
}

#endif // CARTESIAN1DSCALARRMSWRITER_HPP
