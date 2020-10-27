/** 
 * @file ISphericalTorPolMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLMSPECTRUMWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolEnergyBaseWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ISphericalTorPolMSpectrumWriter: public ISphericalTorPolEnergyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolMSpectrumWriter();

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
          * @brief Storage for the Toroidal energy
          */
         Array mTorEnergy;

         /**
          * @brief Storage for the Poloidal energy
          */
         Array mPolEnergy;

         /**
          * @brief Initialize energy storage
          */
         virtual void initializeEnergy();

         /**
          * @brief Store energy from Q component
          */
         virtual void storeQEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from S component
          */
         virtual void storeSEnergy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store energy from T component
          */
         virtual void storeTEnergy(const int l, const int m, const MHDFloat energy);
   };

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLMSPECTRUMWRITER_HPP
