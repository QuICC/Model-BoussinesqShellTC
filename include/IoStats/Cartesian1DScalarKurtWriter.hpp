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
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) skew calculation for a scalar field
    */
   class Cartesian1DScalarKurtWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarKurtWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarKurtWriter& RMS, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarKurtWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         void precompute(Transform::TransformCoordinatorType& coord);
         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         void compute(Transform::TransformCoordinatorType& coord);
         void postcompute(Transform::TransformCoordinatorType& coord);
         
         /**
          * @brief Write State to file
          */
         virtual void prewrite();
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

         /**
          * @brief Storage for the scalar energy
          */
         Array mAvg = Avg;
         Array mRMS = RMS;
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
