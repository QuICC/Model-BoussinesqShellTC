/** 
 * @file CylinderTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII energy calculation for a Toroidal/Poloidal field in a cylinder (Worland expansion)
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef CYLINDERTORPOLENERGYWRITER_HPP
#define CYLINDERTORPOLENERGYWRITER_HPP

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
#include "IoVariable/IVariableAsciiWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII energy calculation for a Toroidal/Poloidal field in a cylinder (Worland expansion)
    */
   class CylinderTorPolEnergyWriter: public IVariableAsciiWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         CylinderTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~CylinderTorPolEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();

         /**
          * @brief Compute energy for toroidal/poloidal field
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
          * @brief Cylindrical volume to normalize energy to energy density
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
          * @brief Worland operator to integrate in radius for energy
          */
         Matrix mRIntgOp;

         /**
          * @brief Chebyshev operator to integrate in Z for energy
          */
         Matrix mZIntgOp;
   };

   inline bool CylinderTorPolEnergyWriter::isHeavy() const
   {
      return true;
   }

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<CylinderTorPolEnergyWriter> SharedCylinderTorPolEnergyWriter;

}
}

#endif // CYLINDERTORPOLENERGYWRITER_HPP
