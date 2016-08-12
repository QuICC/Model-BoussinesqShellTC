/** 
 * @file Cartesian1DScalarAvgWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSCALARAVGWRITER_HPP
#define CARTESIAN1DSCALARAVGWRITER_HPP

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

namespace GeoMHDiSCC {

namespace IoStats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) horizontal average calculation for a scalar field
    */
   class Cartesian1DScalarAvgWriter: public IStatisticsAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarAvgWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarAvgWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Pre computation stage for average of scalar field
          */
         void preCompute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Compute average of scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Post computation stage for average of scalar field
          */
         void postCompute(Transform::TransformCoordinatorType& coord);
         
         /**
          * @brief Write State to file
          */
         virtual void write();

         /**
          * To share average with other stats
          */
         const Array& average() const;

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
         Array mAvg;
   };

   inline bool Cartesian1DScalarAvgWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DScalarAvgWriter> SharedCartesian1DScalarAvgWriter;

}
}
   
#endif // CARTESIAN1DSCALARAVGWRITER_HPP
