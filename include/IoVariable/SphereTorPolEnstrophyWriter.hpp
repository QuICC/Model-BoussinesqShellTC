/** 
 * @file SphereTorPolEnstrophyWriter.hpp
 * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a sphere
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef SPHERETORPOLENSTROPHYWRITER_HPP
#define SPHERETORPOLENSTROPHYWRITER_HPP

// Configuration includes
//
#include "SmartPointers/SharedPtrMacro.h"

// System includes
//

// External includes
//

// Project includes
//
#include "IoVariable/ISphericalTorPolEnstrophyWriter.hpp"

namespace QuICC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII spherical harmonics enstrophy calculation for a Toroidal/Poloidal field in a sphere
    */
   class SphereTorPolEnstrophyWriter: public ISphericalTorPolEnstrophyWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param prefix Prefix to use for file name
          * @param type Type of the file (typically scheme name)
          */
         SphereTorPolEnstrophyWriter(const std::string& prefix, const std::string& type);

         /**
          * @brief Destructor
          */
         virtual ~SphereTorPolEnstrophyWriter();

         /**
          * @brief Initialise the operator, transform and file
          */
         virtual void init();
         
      protected:

      private:
   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<SphereTorPolEnstrophyWriter> SharedSphereTorPolEnstrophyWriter;

}
}

#endif // SPHERETORPOLENSTROPHYWRITER_HPP