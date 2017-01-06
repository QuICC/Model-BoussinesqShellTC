/** 
 * @file ShellScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical shell
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SHELLSCALARENERGYWRITER_HPP
#define SHELLSCALARENERGYWRITER_HPP

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

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a spherical shell
    */
   class ShellScalarEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ShellScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ShellScalarEnergyWriter();

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
         /*
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the scalar energy
          */
         MHDFloat mEnergy;

         /**
          * @brief Chebyshev operator for spherical integral in radius (include r^2 factor)
          */
         SparseMatrix mSphIntgOp;

   };

   inline bool ShellScalarEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<ShellScalarEnergyWriter> SharedShellScalarEnergyWriter;

}
}

#endif // SHELLSCALARENERGYWRITER_HPP
