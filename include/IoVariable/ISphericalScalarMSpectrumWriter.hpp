/** 
 * @file ISphericalScalarMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a scalar field in a spherical geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALSCALARMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALSCALARMSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//
#include <memory>

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalScalarEnergyBaseWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a scalar field in a sphere
    */
   class ISphericalScalarMSpectrumWriter: public ISphericalScalarEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalScalarMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalScalarMSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Storage for the scalar energy
          */
         Array mEnergy;

         /**
          * @brief Initialize energy storage
          */
         virtual void initializeEnergy();

         /**
          * @brief Store energy
          */
         virtual void storeEnergy(const int l, const int m, const MHDFloat energy);
   };

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALSCALARMSPECTRUMWRITER_HPP
