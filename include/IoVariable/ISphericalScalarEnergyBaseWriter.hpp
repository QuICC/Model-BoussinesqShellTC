/** 
 * @file ISphericalScalarEnergyBaseWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYBASEWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYBASEWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "Base/Typedefs.hpp"
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a sphere
    */
   class ISphericalScalarEnergyBaseWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param name       Filename
          * @param ext        File extension
          * @param header     Header string of file
          * @param type       Type string of file
          * @param version    Version string of file
          * @param id         ID of the dimension space
          * @param mode       Write mode of file
          */
         ISphericalScalarEnergyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const WriteMode mode = EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarEnergyBaseWriter();

         /**
          * @brief Activate output of parity splitting in energy output
          */
         void showParity();

         /**
          * @brief Compute energy for scalar field
          */
         void compute(Transform::TransformCoordinatorType& coord);

         /**
          * @brief Requires heavy calculation?
          */
         virtual bool isHeavy() const; 
         
      protected:
         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

         /*
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

         /**
          * @brief Flag to show parity split in energy
          */
         bool mShowParity;

      private:
         /**
          * @brief Initialize energy storage
          */
         virtual void initializeEnergy() = 0;

         /**
          * @brief Store energy
          */
         virtual void storeEnergy(const int l, const int m, const MHDFloat energy) = 0;
   };

   inline bool ISphericalScalarEnergyBaseWriter::isHeavy() const
   {
      return true;
   }

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYBASEWRITER_HPP
