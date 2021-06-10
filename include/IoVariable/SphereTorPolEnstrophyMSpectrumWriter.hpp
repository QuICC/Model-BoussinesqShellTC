/** 
 * @file SphereTorPolEnstrophyMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L enstrophy spectrum calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolEnstrophyMSpectrumWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics M enstrophy spectrum calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnstrophyMSpectrumWriter: public ISphericalTorPolEnstrophyMSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnstrophyMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnstrophyMSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereTorPolEnstrophyMSpectrumWriter> SharedSphereTorPolEnstrophyMSpectrumWriter;

}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLENSTROPHYMSPECTRUMWRITER_HPP
