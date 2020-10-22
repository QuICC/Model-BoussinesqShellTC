/** 
 * @file SphereTorPolMSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLMSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLMSPECTRUMWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolMSpectrumWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics M energy spectrum calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolMSpectrumWriter: public ISphericalTorPolMSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolMSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolMSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereTorPolMSpectrumWriter> SharedSphereTorPolMSpectrumWriter;

}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLMSPECTRUMWRITER_HPP
