/** 
 * @file Cartesian1DScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSCALARENERGYWRITER_HPP
#define CARTESIAN1DSCALARENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a scalar field
    */
   class Cartesian1DScalarEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DScalarEnergyWriter();

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
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

         /**
          * @brief Chebyshev operator for integration
          */
         SparseMatrix mIntgOp;

   };

   inline bool Cartesian1DScalarEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DScalarEnergyWriter> SharedCartesian1DScalarEnergyWriter;

}
}

#endif // CARTESIAN1DSCALARENERGYWRITER_HPP
