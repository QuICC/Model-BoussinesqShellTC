/** 
 * @file SphereTorPolEnergyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERETORPOLENERGYWRITER_HPP
#define SPHERETORPOLENERGYWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolEnergyWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics energy calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnergyWriter: public ISphericalTorPolEnergyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnergyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnergyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereTorPolEnergyWriter> SharedSphereTorPolEnergyWriter;

}
}

#endif // SPHERETORPOLENERGYWRITER_HPP
