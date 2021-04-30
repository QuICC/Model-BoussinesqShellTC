/** 
 * @file ISphericalTorPolEnstrophyBaseWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a spherical geometry
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWBASERITER_HPP
#define QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYWBASERITER_HPP

// Configuration includes
//

// System includes
//

// External includes
//

// Project includes
//
#include "Enums/FieldIds.hpp"
#include "Resolutions/Resolution.hpp"
#include "IoVariable/IVariableAsciiWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a spherical shell
    */
   class ISphericalTorPolEnstrophyBaseWriter: public IVariableAsciiWriter
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
         ISphericalTorPolEnstrophyBaseWriter(std::string name, std::string ext, std::string header, std::string type, std::string version, const Dimensions::Space::Id id, const WriteMode mode = EXTEND);

         /**
          * @brief Destructor
          */
         virtual ~ISphericalTorPolEnstrophyBaseWriter();

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

         /**
          * @brief Spherical volume to normalize enstrophy
          */
         MHDFloat mVolume;

         /**
          * @brief Flag to show parity split in enstrophy
          */
         bool mShowParity;

      private:
         /**
          * @brief Initialize enstrophy storage
          */
         virtual void initializeEnergy() = 0;

         /**
          * @brief Store first toroidal contribution to enstrophy
          */
         virtual void storeT1Enstrophy(const int l, const int m, const MHDFloat energy) = 0;

         /**
          * @brief Store second toroidal contribution to enstrophy
          */
         virtual void storeT2Enstrophy(const int l, const int m, const MHDFloat energy) = 0;

         /**
          * @brief Store poloidal contribution to enstrophy
          */
         virtual void storePEnstrophy(const int l, const int m, const MHDFloat energy) = 0;
   };

   inline bool ISphericalTorPolEnstrophyBaseWriter::isHeavy() const
   {
      return true;
   }

}
}

#endif // QUICC_IO_VARIABLE_ISPHERICALTORPOLENERGYBASEWRITER_HPP
