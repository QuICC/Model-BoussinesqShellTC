/** 
 * @file SphereTorPolLSpectrumWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a Toroidal/Poloidal field in a sphere
 */

#ifndef QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP
#define QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolLSpectrumWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics L energy spectrum calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolLSpectrumWriter: public ISphericalTorPolLSpectrumWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolLSpectrumWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolLSpectrumWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereTorPolLSpectrumWriter> SharedSphereTorPolLSpectrumWriter;

}
}

#endif // QUICC_IO_VARIABLE_SPHERETORPOLLSPECTRUMWRITER_HPP