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
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoStats {

   /**
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) horizontal average calculation for a scalar field
    */
   class Cartesian1DScalarAvgWriter: public IoVariable::IVariableAsciiEWriter
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

         void precompute(Transform::TransformCoordinatorType& coord);
         /**
          * @brief Compute energy for scalar field
          */
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
