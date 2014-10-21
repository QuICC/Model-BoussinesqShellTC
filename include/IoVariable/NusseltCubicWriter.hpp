/** 
 * @file NusseltCubicWriter.hpp
 * @brief Implementation of the ASCII Nusselt number writer for a finite 3D box
 * @author Philippe Marti \<philippe.marti@colorado.edu\>
 */

#ifndef NUSSELTCUBICWRITER_HPP
#define NUSSELTCUBICWRITER_HPP

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
#include "IoVariable/IVariableAsciiEWriter.hpp"
#include "TypeSelectors/ScalarSelector.hpp"

namespace GeoMHDiSCC {

namespace IoVariable {

   /**
    * @brief Implementation of the ASCII Nusselt number writer for a finite 3D box
    */
   class NusseltCubicWriter: public IVariableAsciiEWriter
   {
      public:
         /**
          * @brief Constructor
          *
          * @param type Type of the file (typically scheme name)
          */
         NusseltCubicWriter(std::string type);

         /**
          * @brief Destructor
          */
         virtual ~NusseltCubicWriter();

         /**
          * @brief Write State to file
          */
         virtual void write();
         
      protected:

      private:

   };

   /// Typedef for a shared pointer of a HDF5 state file writer
   typedef SharedPtrMacro<NusseltCubicWriter> SharedNusseltCubicWriter;

}
}

#endif // NUSSELTCUBICWRITER_HPP
