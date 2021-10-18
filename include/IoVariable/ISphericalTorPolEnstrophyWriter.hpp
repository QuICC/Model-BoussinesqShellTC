/** 
 * @file ISphericalTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP

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
         Array mTorEnstrophy;

         /**
          * @brief Storage for the Poloidal enstrophy
          */
         Array mPolEnstrophy;

         /**
          * @brief Initialize enstrophy storage
          */
         virtual void initializeEnstrophy();

         /**
          * @brief Store first toroidal contribution to enstrophy
          */
         virtual void storeT1Enstrophy(const int l, const int m, const MHDFloat enstrophy);

         /**
          * @brief Store second toroidal contribution to enstrophy
          */
         virtual void storeT2Enstrophy(const int l, const int m, const MHDFloat enstrophy);

         /**
          * @brief Store poloidal contribution to enstrophy
          */
         virtual void storePEnstrophy(const int l, const int m, const MHDFloat enstrophy);
   };

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENSTROPHYWRITER_HPP