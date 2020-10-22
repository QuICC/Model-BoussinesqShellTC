/** 
 * @file ISphericalScalarEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a scalar field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP

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
   class ISphericalScalarEnergyWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarEnergyWriter();

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
         /**
          * @brief Data ordering is m slowest
          */
         bool mHasMOrdering;

         /*
          * @brief Spherical volume to normalize energy to energy density
          */
         MHDFloat mVolume;

      private:

         /**
          * @brief Storage for the scalar energy
          */
         MHDFloat mEnergy;
   };

   inline bool ISphericalScalarEnergyWriter::isHeavy() const
   {
      return true;
   }

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARENERGYWRITER_HPP
