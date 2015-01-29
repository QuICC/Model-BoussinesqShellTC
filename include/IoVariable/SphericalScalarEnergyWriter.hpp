/** 
 * @file SphericalScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERICALSCALARENERGYWRITER_HPP
#define SPHERICALSCALARENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field
    */
   class SphericalScalarEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphericalScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphericalScalarEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

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
         /**
          * @brief Storage for the scalar energy
          */
         MHDFloat mEnergy;

         /**
          * @brief Chebyshev operator to integrate in radius
          */
         SparseMatrix mIntgOp;

   };

   inline bool SphericalScalarEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphericalScalarEnergyWriter> SharedSphericalScalarEnergyWriter;

}
}

#endif // SPHERICALSCALARENERGYWRITER_HPP
