/** 
 * @file Cartesian1DPrimitiveEnergyWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (primitive formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DPRIMITIVEENERGYWRITER_HPP
#define CARTESIAN1DPRIMITIVEENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII Cartesian 1D (double periodic)energy calculation for a vector field (primitive formulation)
    */
   class Cartesian1DPrimitiveEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DPrimitiveEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DPrimitiveEnergyWriter();

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
          * @brief Cartesian box volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the X component energy
          */
         MHDFloat mXEnergy;

         /**
          * @brief Storage for the Y component energy
          */
         MHDFloat mYEnergy;

         /**
          * @brief Storage for the Z component energy
          */
         MHDFloat mZEnergy;

         /**
          * @brief Chebyshev operator to integrate
          */
         SparseMatrix mIntgOp;
   };

   inline bool Cartesian1DPrimitiveEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DPrimitiveEnergyWriter> SharedCartesian1DPrimitiveEnergyWriter;

}
}

#endif // CARTESIAN1DPRIMITIVEENERGYWRITER_HPP
