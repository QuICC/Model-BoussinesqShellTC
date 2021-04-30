/** 
 * @file ISphericalTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Resolutions/Resolution.hpp"
#include "IoVariable/ISphericalTorPolEnstrophyBaseWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a sphere
    */
   class ISphericalTorPolEnstrophyWriter: public ISphericalTorPolEnstrophyBaseWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         ISphericalTorPolEnstrophyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnstrophyWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:
         /**
          * @brief Storage for the Toroidal enstrophy
          */
         Array mTorEnergy;

         /**
          * @brief Storage for the Poloidal enstrophy
          */
         Array mPolEnergy;

         /**
          * @brief Initialize enstrophy storage
          */
         virtual void initializeEnergy();

         /**
          * @brief Store first toroidal contribution to enstrophy
          */
         virtual void storeT1Enstrophy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store second toroidal contribution to enstrophy
          */
         virtual void storeT2Enstrophy(const int l, const int m, const MHDFloat energy);

         /**
          * @brief Store poloidal contribution to enstrophy
          */
         virtual void storePEnstrophy(const int l, const int m, const MHDFloat energy);
   };

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWRITER_HPP
