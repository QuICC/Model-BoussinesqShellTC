/** 
 * @file SphericalTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERICALTORPOLENERGYWRITER_HPP
#define SPHERICALTORPOLENERGYWRITER_HPP

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
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field
    */
   class SphericalTorPolEnergyWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphericalTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphericalTorPolEnergyWriter();

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

   inline bool SphericalTorPolEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphericalTorPolEnergyWriter> SharedSphericalTorPolEnergyWriter;

}
}

#endif // SPHERICALTORPOLENERGYWRITER_HPP
