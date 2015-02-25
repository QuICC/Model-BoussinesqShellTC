/** 
 * @file Cartesian1DStreamEnergyWriter.hpp
 * @brief Implementation of the ASCII Cartesian 1D (double periodic) energy calculation for a vector field (streamfunction formulation)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CARTESIAN1DSTREAMENERGYWRITER_HPP
#define CARTESIAN1DSTREAMENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII Cartesian 1D (double periodic)energy calculation for a vector field (streamfunction formulation)
    */
   class Cartesian1DStreamEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         Cartesian1DStreamEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~Cartesian1DStreamEnergyWriter();

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
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Storage for the Toroidal energy
          */
         MHDFloat mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         MHDFloat mPolEnergy;

         /**
          * @brief Chebyshev operator to integrate in radius
          */
         SparseMatrix mIntgOp;

         /**
          * @brief Chebyshev operator for spherical integral in radius (include r^2 factor)
          */
         SparseMatrix mSphIntgOp;
   };

   inline bool Cartesian1DStreamEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<Cartesian1DStreamEnergyWriter> SharedCartesian1DStreamEnergyWriter;

}
}

#endif // CARTESIAN1DSTREAMENERGYWRITER_HPP
