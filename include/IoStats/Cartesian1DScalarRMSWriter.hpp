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
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
    */
   class Cartesian1DScalarRMSWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarRMSWriter(const std::string& prefix, const SharedCartesian1DScalarAvgWriter& Avg, const SharedCartesian1DScalarRMSWriter& RMS, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarRMSWriter();

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
         /*
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mArea;

         /**
          * @brief Storage for the scalar energy
          */
         Array mAvg = Avg;
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
